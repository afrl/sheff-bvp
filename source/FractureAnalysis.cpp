#include "fea_analysis/FractureAnalysis.h"

#include "flex_mesh/maps/AdjacencyLists.h"
#include "flex_mesh/algorithms/Search.h"

#include "simple_math/Utility.h"

#include "finite_element/FiniteElement.h"
#include "finite_element/BaseFormulation.h"
#include "finite_element/elasticity/InfinitesimalStrain.h"
#include "finite_element/damage/XFEMFormulation.h"

#include "fea_crack/FeaCrackUtils.h"

#include "sparse_linear_algebra/Factory.h"

#include <memory>

namespace fea_analysis
{
    FractureAnalysis::FractureAnalysis()
    {
    }

    void FractureAnalysis::setup_analysis()
    {
        LinearElasticAnalysis::setup_analysis();
        adjacency_list_.initialize(*mesh_);
        assembler_->set_element_workspaces();
        enrichment_surfaces_ = std::make_unique<fea_crack::EnrichmentSurfaceManager>(
            fea_crack::RepresentationType::LevelSet,
            *dof_map_, *mesh_, adjacency_list_);
    }

    void FractureAnalysis::insert_precracks(const std::vector<simple_math::Point>& initiation_points)
    {
        // Initiate all of the cracks
        for (auto i : utility::indices(initiation_points))
        {
            auto& new_crack = enrichment_surfaces_->add_enrichment_surface();
            auto initial_element_ptr = flex_mesh::find_element_with_point(*mesh_, initiation_points[i]);
            if (!initial_element_ptr)
            {
                throw _UTIL_MSG_EXCEPT("No element was found that contains the point "
                    + simple_math::to_string(initiation_points[i]));
            }
            auto& direction_criteria = initial_element_ptr->get_direction_criterion();
            simple_math::VoigtVector zero_stress = simple_math::VoigtVector::Zero(6);
            auto crack_normal = direction_criteria.eval(i, initiation_points[i], zero_stress);
            new_crack.initiate(initial_element_ptr, 0, initiation_points[i], crack_normal);

            //cracks.back()->write_vtk("./", 0);
        }

        // Grow all the cracks one element until they can no longer grow
        int growth_step = 1;
        auto element_num_dofs = dof_map_->get_array_of_element_n_dofs();
        for (auto& crack : *enrichment_surfaces_)
        {
            auto& crack_state = crack.get_crack_state();
            while (!crack_state.get_enrichments_at_front().empty())
            {
                std::vector<fea_crack::EnrichmentSurfaceState::EnrichmentAtFront*> enrichments_at_front_to_grow;
                for (auto& enrichment_at_front : crack_state.get_enrichments_at_front())
                {
                    auto element_ptr = enrichment_at_front->get_enrichment()->get_element();
                    auto& criterion = element_ptr->get_propagation_criterion();
                    if (criterion.eval(crack_state.get_enrichment_surface_id(),
                        element_ptr->get_mean_coordinate(), simple_math::Vector3D::Zero()))
                    {
                        enrichments_at_front_to_grow.push_back(enrichment_at_front.get());
                    }
                }

                // If no enrichments can grow due to the propagation criteria, then stop
                if (enrichments_at_front_to_grow.empty())
                {
                    break;
                }

                for (auto& enrichment_at_front_to_grow : enrichments_at_front_to_grow)
                {
                    std::vector<fea_crack::LeafDofInformation> leaves_to_grow_into;
                    std::vector<simple_math::GeneralVector> solution_of_ahead_of_tip_elements;
                    auto& elements_ahead_of_tip = enrichment_at_front_to_grow->get_adjacent_elements_ahead_of_front();
                    for (auto& element_ahead_of_tip : elements_ahead_of_tip)
                    {
                        auto ahead_of_tip_element_index = element_ahead_of_tip->get_element()->get_local_number();
                        solution_of_ahead_of_tip_elements.emplace_back(element_num_dofs[ahead_of_tip_element_index]);
                        solution_of_ahead_of_tip_elements.back().setZero();
                    }

                    crack.grow_crack_front(*enrichment_at_front_to_grow,
                        solution_of_ahead_of_tip_elements);
                }
                //crack->write_vtk("./", growth_step);
                growth_step++;
            }
        }

        for (auto& crack : *enrichment_surfaces_)
        {
            crack.write_vtk("./", 0, 0);
        }
    }

    void FractureAnalysis::solve_time_step()
    {
        gl_stiffness_->update_local_global_dof_map(
            mesh_->elements().size(),
            dof_map_->get_n_local_dofs(),
            dof_map_->get_n_global_dofs());
        gl_displacement_->update_local_global_dof_map(
            mesh_->elements().size(),
            dof_map_->get_n_local_dofs(),
            dof_map_->get_n_global_dofs());
        gl_eigenstrain_force_->update_local_global_dof_map(
            mesh_->elements().size(),
            dof_map_->get_n_local_dofs(),
            dof_map_->get_n_global_dofs());
        gl_rhs_->update_local_global_dof_map(
            mesh_->elements().size(),
            dof_map_->get_n_local_dofs(),
            dof_map_->get_n_global_dofs());

        apply_boundary_conditions();

        solve_equillibrium_step();

        /*
        // TODO indicate if Secant or Tangent with BCUpdateType::Secant
        apply_boundary_conditions();

        solve_equillibrium_step();

        bool found_equillibrium = false;

        while (!found_equillibrium)
        {
            solve_equillibrium_step();

            auto growth_successful = try_to_grow_cracks();

            if (!growth_successful)
            {
                auto initations_successful = try_to_initiate_cracks();
                found_equillibrium = !initations_successful;
            }
            else
            {
                found_equillibrium = false;
            }
        }
        */
    }

    void FractureAnalysis::solve_equillibrium_step()
    {
        //////////////
        // Assembly //
        //////////////

        // Queue assembly of the global stiffness matrix
        assembler_->queue_matrix_assembly(
            gl_stiffness_.get(),
            [](finite_element::FiniteElement & fe) -> auto
        {
            auto element_matrix = fe.get_formulation().get_element_matrix();

            std::string output_path = "ElementStiffness_" + std::to_string(fe.get_global_number()) + ".txt";
            std::ofstream output_stream(output_path);
            output_stream << element_matrix.str() << std::endl;
            output_stream.close();

            return element_matrix;
        });

        {
            std::string output_path = "DoFMaps.txt";
            std::ofstream output_stream(output_path);
            output_stream << dof_map_->str() << std::endl;
            output_stream.close();
        }

        // Queue assembly of the global vector containing RHS terms pertaining
        // to eigenstrains (such as thermal strain, plastic strain, etc.)
        assembler_->queue_vector_assembly(
            gl_eigenstrain_force_.get(),
            [](finite_element::FiniteElement& fe) -> auto
        {
            return fe.get_formulation().get_element_rhs_vector();
        });

        // Execute both assemblies as one element loop
        assembler_->execute(gl_displacement_.get());

        // Total RHS = constraint RHS + applied loads + inelastic RHS (temperature/moisture terms)
        *gl_rhs_ = *gl_rhs_ + *gl_eigenstrain_force_;

        ///////////
        // Solve //
        ///////////

        // TODO make the solver automatically update the ghost values of the solution
        sparse_solver_->solve(
            gl_stiffness_,
            gl_rhs_,
            gl_displacement_);

        // Debug
        auto solution = gl_displacement_->get_as_eigen_vector();
        fea_crack::write_enrichment_surfaces_vtk_tmp(
            *mesh_, *dof_map_, "Solution.vtk", solution);
    }

    bool FractureAnalysis::try_to_initiate_cracks()
    {
        bool successful = false;
        std::unique_ptr<finite_element::ElementWorkspace> element_workspace_ptr
            = std::make_unique<finite_element::ElementWorkspace>();
        for (auto& e_ptr : mesh_->elements())
        {
            // Provide the current solution in case the element requries it.
            auto elem_solution_vector = gl_displacement_->get_element_part(e_ptr->get_local_number());

            // This should be separated out in the future, but for now provide the workspace
            e_ptr->set_workspace(element_workspace_ptr.get());

            // Start calculation block
            e_ptr->start_calculations(elem_solution_vector);

            auto& formulation = e_ptr->get_formulation();
            auto elastic_formulation = dynamic_cast<finite_element::InfinitesimalStrainFormulation*>(&formulation);
            if (elastic_formulation)
            {
                auto initiation_result = elastic_formulation->eval_initiation_criterion();
                if (initiation_result.criterion_satisified)
                {
                    auto& new_crack = enrichment_surfaces_->add_enrichment_surface();
                    new_crack.initiate(e_ptr.get(), 0, initiation_result.point, initiation_result.crack_plane_normal);
                    successful = true;
                }
            }

            auto rxfem_formulation = dynamic_cast<finite_element::XFEMFormulation*>(&formulation);
            if (rxfem_formulation)
            {
                auto initiation_results = rxfem_formulation->eval_initiation_criterion();
                for (auto i : utility::indices(initiation_results))
                {
                    auto& initiation_result = initiation_results[i];
                    if (initiation_result.criterion_satisified)
                    {
                        auto& new_crack = enrichment_surfaces_->add_enrichment_surface();
                        new_crack.initiate(e_ptr.get(), i, initiation_result.point, initiation_result.crack_plane_normal);
                        successful = true;
                    }
                }
            }

            // End calculation block
            e_ptr->end_calculations();
        }
        return successful;
    }

    bool FractureAnalysis::try_to_grow_cracks()
    {
        /*
        bool growth_successful = false;
        for (auto& crack : cracks_)
        {
            auto& crack_state = crack->get_crack_state();
            auto& enrichments_at_front = crack_state.get_enrichments_at_front();

            std::vector<fea_crack::CrackState::Index> enrichments_at_crack_tip_to_grow;
            if (enrichments_at_front.size() > 0)
            {
                enrichments_at_crack_tip_to_grow.push_back(0);
            }

            auto element_num_dofs = dof_map_->get_array_of_element_n_dofs();
            for (auto& enrichment_to_grow_index : enrichments_at_crack_tip_to_grow)
            {
                std::vector<fea_crack::LeafDofInformation> leaves_to_grow_into;
                std::vector<simple_math::GeneralVector> solution_of_ahead_of_tip_elements;
                auto& enrichment_to_grow = crack_state.get_enrichments_at_front()[enrichment_to_grow_index];
                auto& elements_ahead_of_tip = enrichment_to_grow->get_adjacent_elements_ahead_of_front();
                for (auto& element_ahead_of_tip : elements_ahead_of_tip)
                {
                    auto ahead_of_tip_element_index = element_ahead_of_tip->get_element()->get_local_number();
                    solution_of_ahead_of_tip_elements.emplace_back(element_num_dofs[ahead_of_tip_element_index]);
                    solution_of_ahead_of_tip_elements.back().setZero();
                }

                crack->grow_crack_front(enrichment_to_grow_index, solution_of_ahead_of_tip_elements);
                growth_successful = true;
            }
        }

        return growth_successful;
        */
        return false;
    }
}