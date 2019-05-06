#include "fea_analysis/FractureAnalysis.h"

#include "flex_mesh/maps/AdjacencyLists.h"

#include "simple_math/Utility.h"

#include "finite_element/FiniteElement.h"
#include "finite_element/BaseFormulation.h"
#include "finite_element/elasticity/InfinitesimalStrain.h"

#include "fea_crack/RXFEMCrack.h"

#include "sparse_linear_algebra/Factory.h"

namespace fea_analysis
{
    void FractureAnalysis::solve_time_step()
    {
        // TODO indicate if Secant or Tangent with BCUpdateType::Secant
        apply_boundary_conditions();

        // As a test, initiate a crack in element 0
        {
            fea_crack::RXFEMCrack::ElementAdjacencyList adjacency_list;
            adjacency_list.initialize(*mesh_);

            fea_crack::RXFEMCrack crack(0, fea_crack::CrackRepType::LevelSet,
                *dof_map_, *mesh_, adjacency_list);

            simple_math::Point initial_coord = { 0.5, 0.5, 0.5 };
            simple_math::Point initial_norm = { 1.0, 0.0, 0.0 };
            simple_math::Vector3D initial_normal = simple_math::point_as_vector(initial_norm);

            // Initiate the crack
            finite_element::ElementWorkspace e_workspace;
            mesh_->elements().front()->set_workspace(&e_workspace);
            crack.initiate(mesh_->elements().front(), 0, initial_coord, initial_normal);

            // Call update_local_global_dof_map for each container
            const auto n_elements = static_cast<flex_mesh::LocalIndex>(mesh_->elements().size());
            gl_displacement_->update_local_global_dof_map(n_elements,
                dof_map_->get_n_local_dofs(), dof_map_->get_n_global_dofs());
            gl_eigenstrain_force_->update_local_global_dof_map(n_elements,
                dof_map_->get_n_local_dofs(), dof_map_->get_n_global_dofs());
            gl_rhs_->update_local_global_dof_map(n_elements,
                dof_map_->get_n_local_dofs(), dof_map_->get_n_global_dofs());
            gl_stiffness_->update_local_global_dof_map(n_elements,
                dof_map_->get_n_local_dofs(), dof_map_->get_n_global_dofs());
        }

        //////////////
        // Assembly //
        //////////////

        // Queue assembly of the global stiffness matrix
        assembler_->queue_matrix_assembly(
            gl_stiffness_.get(),
            [](finite_element::FiniteElement& fe) -> auto
        {
            return fe.get_formulation().get_element_matrix();
        });

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
    }
}