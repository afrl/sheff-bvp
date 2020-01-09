#include "fea_analysis/FractureAnalysis.h" 
#include "flex_mesh/maps/AdjacencyLists.h"
#include "flex_mesh/algorithms/Search.h"
#include "flex_mesh/element/PolyhedronElement.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/IOFeaIo.h"

#include "simple_math/Utility.h"

#include "finite_element/FiniteElement.h"
#include "finite_element/BaseFormulation.h"
#include "finite_element/elasticity/InfinitesimalStrain.h"
#include "finite_element/damage/XFEMFormulation.h"

#include "sparse_linear_algebra/LinAlgFactory.h"

#include "output/Output.h"
#include "output/OutputFeaAnalysisHelpers.h"

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

        // Create the callback for the fea_crack module to get the element solution
        // vector when evaluating the direction criteria.
        fea_crack::GetElementSolutionCallback get_element_solution_cb =
        [gl_displacement = gl_displacement_.get()](const Element& element) -> simple_math::GeneralVector
        {
            return gl_displacement->get_element_part(element.get_local_number());
        };

        enrichment_surfaces_ = std::make_unique<fea_crack::EnrichmentSurfaceManager>(
            *mesh_,
            adjacency_list_,
            get_element_solution_cb,
            *dof_map_);

    }

    void FractureAnalysis::insert_precracks(const std::vector<simple_math::Point>& initiation_points)
    {
        dof_map_->start_enrichment_phase();

        // Initiate all of the cracks
        for (auto i : utility::indices(initiation_points))
        {
            auto new_crack_ptr = enrichment_surfaces_->add_enrichment_surface();
            auto initial_element_ptr = flex_mesh::find_element_with_point(*mesh_, initiation_points[i]);
            if (!initial_element_ptr)
            {
                throw _UTIL_MSG_EXCEPT("No element was found that contains the point "
                    + simple_math::to_string(initiation_points[i]));
            }
            auto& direction_criteria = initial_element_ptr->get_direction_criterion();
            simple_math::VoigtVector zero_stress = simple_math::VoigtVector::Zero(6);
            auto crack_normal = direction_criteria.eval(new_crack_ptr->get_global_id(), initiation_points[i], zero_stress);
            new_crack_ptr->initiate(*initial_element_ptr, 0, initiation_points[i], crack_normal);

            new_crack_ptr->write_vtk("./", 0, 0);
        }

        // Grow all the cracks one element until they can no longer grow
        int growth_step = 1;
        for (auto& crack : *enrichment_surfaces_)
        {
            auto& crack_state = crack.get_representation().get_state();
            auto& enrichments_at_front = crack.get_enrichments_at_front();
            while (!enrichments_at_front.empty())
            {
                std::vector<fea_crack::EnrichmentAtFront*> enrichments_at_front_to_grow;
                for (auto& enrichment_at_front : enrichments_at_front)
                {
                    auto& element = enrichment_at_front->get_enrichment_ref().get_element();
                    auto& criterion = element.get_propagation_criterion();
                    if (criterion.eval(crack_state.get_enrichment_surface_id(),
                        element.get_mean_coordinate(), simple_math::Vector3D::Zero()))
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
                    auto search_result = crack.search_ahead_of_front(
                        *enrichment_at_front_to_grow);

                    crack.grow_crack_front_into_all_leaves(
                        *enrichment_at_front_to_grow,
                        search_result);
                }
                crack.write_vtk("./", 0, growth_step);
                growth_step++;
            }
        }

        for (auto& crack : *enrichment_surfaces_)
        {
            crack.write_vtk("./", 0, 0);
        }

        // Change phase to static phase
        dof_map_->start_renumber_phase();

        auto& renumbering_new_rank_to_old_rank_map =
            dof_map_->get_renumbering_new_rank_to_old_rank_map();
        auto& renumbering_old_rank_to_global_map =
            dof_map_->get_renumbering_old_rank_to_global_map();

        gl_stiffness_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_displacement_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_rhs_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_eigenstrain_force_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);

        dof_map_->end_renumber_phase();
    }

    void FractureAnalysis::insert_precracks_by_function(
        const std::vector<LevelSetFunction>& level_set_functions)
    {
        dof_map_->start_enrichment_phase();

        for (auto& level_set_function : level_set_functions)
        {
            auto new_surface_ptr = enrichment_surfaces_->add_enrichment_surface();
            new_surface_ptr->insert_based_on_level_set(
                level_set_function,
                [](const simple_math::Point&) { return true; });
        }

        for (auto& crack : *enrichment_surfaces_)
        {
            crack.write_vtk("./", 0, 0);
        }

        // Change phase to static phase
        dof_map_->start_renumber_phase();

        auto& renumbering_new_rank_to_old_rank_map =
            dof_map_->get_renumbering_new_rank_to_old_rank_map();
        auto& renumbering_old_rank_to_global_map =
            dof_map_->get_renumbering_old_rank_to_global_map();

        gl_stiffness_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_displacement_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_rhs_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_eigenstrain_force_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);

        dof_map_->end_renumber_phase();
    }

    void FractureAnalysis::insert_precracks_by_function_outside_previous(
        const std::vector<LevelSetFunction>& level_set_functions)
    {
        dof_map_->start_enrichment_phase();

        for (auto i : utility::indices(level_set_functions))
        {
            auto new_surface_ptr = enrichment_surfaces_->add_enrichment_surface();
            new_surface_ptr->insert_based_on_level_set(
                level_set_functions[i],
                [&level_set_functions, &i](const simple_math::Point& p)
                {
                    for (auto j : utility::range_from_zero(i))
                    {
                        if (level_set_functions[j](p) <= 0.0)
                        {
                            return false;
                        }
                    }
                    return true;
                },
                false);
        }

        for (auto& crack : *enrichment_surfaces_)
        {
            crack.write_vtk("./", 0, 0);
        }

        // Change phase to static phase
        dof_map_->start_renumber_phase();

        auto& renumbering_new_rank_to_old_rank_map =
            dof_map_->get_renumbering_new_rank_to_old_rank_map();
        auto& renumbering_old_rank_to_global_map =
            dof_map_->get_renumbering_old_rank_to_global_map();

        gl_stiffness_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_displacement_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_rhs_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_eigenstrain_force_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);

        dof_map_->end_renumber_phase();
    }

    void FractureAnalysis::insert_precracks_by_function_with_boolean(
        const std::vector<LevelSetFunction>& level_set_functions)
    {
        // TODO relax this with further development.
        if (level_set_functions.size() != 2)
        {
            // Throw error
        }

        // TODO
        auto& f1 = level_set_functions[0];
        auto& f2 = level_set_functions[1];
        LevelSetFunction f1_wrapped =
            [&f1, &f2](const simple_math::Point& p) -> const double
            {
                const double f2_val = f2(p);
                if (f2_val > 0)
                {
                    return f1(p);
                }
                else
                {
                    return (f1(p) - f2_val);
                }
            };

        dof_map_->start_enrichment_phase();

        // Insert surface for f1 (including boolean)
        {
            auto new_surface_ptr = enrichment_surfaces_->add_enrichment_surface();
            new_surface_ptr->insert_based_on_level_set(
                f1_wrapped,
                [](const simple_math::Point& p) { return true; });
        }

        // Insert surface for f2
        {
            // Why do I have to flipped the second condition?  Get to the bottom of this.
            LevelSetFunction f1_wrapped2 =
                [&f1, &f2](const simple_math::Point& p) -> const double
            {
                const double f2_val = f2(p);
                if (f2_val > 0)
                {
                    return f1(p);
                }
                else
                {
                    return (f2_val - f1(p));
                }
            };

            auto new_surface_ptr = enrichment_surfaces_->add_enrichment_surface();
            new_surface_ptr->insert_based_on_level_set(
                f2,
                [f1, f1_wrapped](const simple_math::Point& p) -> const bool {
                    const double f1_val = f1(p);
                    const double f1_wrapped_val = f1_wrapped(p);
                    if (f1_val > 0)
                    {
                        return true;
                    }
                    else if (f1_wrapped_val > 0)
                    {
                        return true;
                    }
                    return false;
                },
                false);
        }

        for (auto& crack : *enrichment_surfaces_)
        {
            crack.write_vtk("./", 0, 0);
        }

        // Change phase to static phase
        dof_map_->start_renumber_phase();

        auto& renumbering_new_rank_to_old_rank_map =
            dof_map_->get_renumbering_new_rank_to_old_rank_map();
        auto& renumbering_old_rank_to_global_map =
            dof_map_->get_renumbering_old_rank_to_global_map();

        gl_stiffness_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_displacement_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_rhs_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);
        gl_eigenstrain_force_->renumber_dofs(0,
            renumbering_old_rank_to_global_map,
            renumbering_new_rank_to_old_rank_map);

        dof_map_->end_renumber_phase();
    }

    void FractureAnalysis::solve_time_step()
    {
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

    static void extract_enriched_elements(
        const Mesh& original_mesh,
        flex_mesh::MeshPartition<Node, flex_mesh::Element<Node>>& enriched_elements_mesh,
        flex_mesh::MeshPartition<Node, flex_mesh::Element<Node>>& enrichment_surfaces_mesh)
    {
        flex_mesh::LocalIndex curr_local_node_i = 0;
        flex_mesh::LocalIndex curr_local_element_i = 0;
        flex_mesh::LocalIndex curr_surf_local_node_i = 0;
        flex_mesh::LocalIndex curr_surf_local_element_i = 0;
        for (auto& e_ptr : original_mesh.elements())
        {
            if (e_ptr->has_formulation_type<finite_element::XFEMFormulation>())
            {
                auto& xfem_form =
                    static_cast<const finite_element::XFEMFormulation&>(
                        e_ptr->get_formulation_const());
                auto tree_traverser = xfem_form.get_tree_traverser();
                auto& integration_domains = xfem_form.get_integration_domains();
                auto& points = integration_domains.get_points();

                flex_mesh::NodeX3Matrix node_coords;
                e_ptr->get_node_coordinates(node_coords);

                // Add nodes and elements for the integration domain of each field
                // Fills enriched_elements_mesh
                for (auto leaf_index : utility::range_from_zero(xfem_form.get_n_leaves()))
                {
                    auto leaf_info = tree_traverser.get_leaf(leaf_index);
                    auto& polyhedron = integration_domains.get_polyhedra()[leaf_info.field_index];
                    auto point_indices = polyhedron.get_point_indices();

                    // Add points for polyhedra
                    flex_mesh::Connectivity<Node, Element::ConnSize> nodes_for_points;
                    std::map<Element::ConnSize, Element::ConnSize> poly_to_facet_point_index_map;
                    Element::ConnSize curr_facet_index = 0;
                    for (auto i : point_indices)
                    {
                        // The points of the integration domain is in the parametric coordinate
                        // system, so they need to be transformed to the global for plotting.
                        auto p = integration_domains.get_points()[i];
                        auto shape_function_values =
                            e_ptr->get_shape_function_values(p);

                        simple_math::Point global_p;
                        for (auto i = 0; i < 3; i++)
                        {
                            global_p[i] = shape_function_values.dot(node_coords.col(i));
                        }

                        nodes_for_points.push_back(
                            enriched_elements_mesh.add_node_by_coord(global_p));
                        nodes_for_points.back()->set_local_number(curr_local_node_i++);
                        poly_to_facet_point_index_map[i] = curr_facet_index;
                        curr_facet_index++;
                    }

                    std::vector<utility::DynamicArray<Element::ConnSize, Element::ConnSize>> indices_for_facets;
                    for (auto& facet : polyhedron.get_facets())
                    {
                        indices_for_facets.emplace_back();
                        indices_for_facets.back().resize(facet.point_indices.size());
                        auto curr_facet_i = 0;
                        for (auto i : facet.point_indices)
                        {
                            indices_for_facets.back()[curr_facet_i] = poly_to_facet_point_index_map[i];
                            curr_facet_i++;
                        }
                    }

                    auto poly_element = std::make_unique<flex_mesh::PolyhedronElement<Node>>();
                    poly_element->initialize(nodes_for_points, indices_for_facets);
                    poly_element->set_local_number(curr_local_element_i++);
                    poly_element->set_global_number(e_ptr->get_global_number());
                    poly_element->set_component_id(e_ptr->get_component_id());
;                   enriched_elements_mesh.elements().add(std::move(poly_element));
                }

                // Add nodes and elements for each enrichment
                // Fills enrichment_surfaces_mesh
                auto& enrichment_facets = integration_domains.get_enrichment_facets();
                for (auto enrichment_index : utility::range_from_zero(xfem_form.get_n_enrichments()))
                {
                    auto& facet_point_indices = enrichment_facets[enrichment_index];
                    
                    std::vector<Node*> connectivity;
                    for (auto i : facet_point_indices)
                    {
                        // The points of the integration domain is in the parametric coordinate
                        // system, so they need to be transformed to the global for plotting.
                        auto p = integration_domains.get_points()[i];
                        auto shape_function_values =
                            e_ptr->get_shape_function_values(p);

                        simple_math::Point global_p;
                        for (auto i = 0; i < 3; i++)
                        {
                            global_p[i] = shape_function_values.dot(node_coords.col(i));
                        }

                        connectivity.push_back(
                            enrichment_surfaces_mesh.add_node_by_coord(global_p));
                        connectivity.back()->set_local_number(curr_surf_local_node_i++);
                    }

                    auto polygon_element = std::make_unique<flex_mesh::Element<Node>>();
                    polygon_element->initialize(
                        flex_mesh::ElementType::POLYGON,
                        connectivity,
                        flex_mesh::LCSAssignment::None);

                    polygon_element->set_local_number(curr_surf_local_element_i++);
                    polygon_element->set_global_number(e_ptr->get_global_number());
                    polygon_element->set_component_id(e_ptr->get_component_id());
                    // Right now, this throws off the reader since the component selection
                    // does not apply to enrichment surfaces.  How do we want to indicate surface?
                        //tree_traverser.get_enrichment(enrichment_index).enrichment_surface_id);
                    enrichment_surfaces_mesh.elements().add(std::move(polygon_element));
                }
            }
        }
    }

    void FractureAnalysis::write_output()
    {
        LinearElasticAnalysis::write_output();

        auto& output_meta_info = output_meta_.get_analysis_meta();
        auto& time_step_info = output_meta_.get_time_step(
            time_step_manager_->get_time_step(),
            time_step_manager_->get_current_time());

        if (output_meta_info.levels_of_detail_names.size() == 1)
        {
            output_meta_info.levels_of_detail_names.push_back("enriched");
        }

        // Update the output meta info with the mesh object for the enriched elements
        fs::path enriched_elements_path_stem = curr_time_step_dir_ / "enriched_elements";
        fs::path enrichment_surfaces_path_stem = curr_time_step_dir_ / "enrichment_surfaces";

        // Note: there is an assumption here that the time step directory is the same folder
        // as the output meta file.  This may need to be generalized.
        std::string enriched_elements_rel_stem = 
            (curr_time_step_dir_.filename() / "enriched_elements").generic_string();
        std::string enrichment_surfaces_rel_stem =
            (curr_time_step_dir_.filename() / "enrichment_surfaces").generic_string();

        fea_io::MeshObject enriched_elements_mesh_object;
        enriched_elements_mesh_object.type = fea_io::MeshObjectType::Volume;
        enriched_elements_mesh_object.mesh_path = enriched_elements_rel_stem;
		enriched_elements_mesh_object.data_paths.push_back(enriched_elements_rel_stem);

        fea_io::MeshObject enrichment_surfaces_mesh_object;
        enrichment_surfaces_mesh_object.type = fea_io::MeshObjectType::Surface;
        enrichment_surfaces_mesh_object.mesh_path = enrichment_surfaces_rel_stem;

        // Update the time step info
        time_step_info.levels_of_detail.push_back({
            1, // level
            fea_io::LODParentRelation::ConformingOverlapping, // relationship_with_parent
            { enriched_elements_mesh_object,
              enrichment_surfaces_mesh_object } }); // mesh_objects

        {
            // Create a mesh for the enriched elements and enrichment surfaces which is only topological
            flex_mesh::MeshPartition<Node, flex_mesh::Element<Node>> enriched_elements_mesh;
            flex_mesh::MeshPartition<Node, flex_mesh::Element<Node>> enrichment_surfaces_mesh;
            extract_enriched_elements(*mesh_,
                enriched_elements_mesh,
                enrichment_surfaces_mesh);

            // DEBUG TODO delete
            flex_mesh::write_vtk(enriched_elements_mesh,
                enriched_elements_path_stem.string() + ".vtk", false, false);
            flex_mesh::write_vtk(enrichment_surfaces_mesh,
                enrichment_surfaces_path_stem.string() + ".vtk", false, false);

            // Write the mesh file that the updated meta info refers to
            flex_mesh::FEAIOInterface<flex_mesh::MeshPartition<Node, flex_mesh::Element<Node>>>::write(
                enriched_elements_mesh, enriched_elements_path_stem.string(), {});
            flex_mesh::FEAIOInterface<flex_mesh::MeshPartition<Node, flex_mesh::Element<Node>>>::write(
                enrichment_surfaces_mesh, enrichment_surfaces_path_stem.string(), {});

            // Write data for enriched elements
            // TODO: write data for enrichment surfaces
            {
                // Create a writer for this time step
                fea_io::DataWriter<flex_mesh::GlobalIndex, flex_mesh::LocalIndex, flex_mesh::FloatType> writer;

                std::vector<output::OutputFlag*> flags_to_be_output;
                for (auto& flag : output_flags_)
                {
                    if (flag->should_be_output())
                    {
                        flags_to_be_output.push_back(flag.get());
                    }
                }

                // TODO only update this if the mesh has changed in subsequent writes
                flex_mesh::SubcomponentMaps enriched_elements_subcomponent_maps;
                enriched_elements_subcomponent_maps.update_maps(enriched_elements_mesh);

                std::ofstream ascii_writer((enriched_elements_path_stem / ".txt").string());

                // Create an output algorithm for each output flag involved in this write
                auto output_algorithms = output::get_output_algorithms(ascii_writer, writer,
                    enriched_elements_mesh.get_partition_id(), flags_to_be_output,
                    static_cast<flex_mesh::LocalIndex>(enriched_elements_mesh.nodes().size()),
                    subcomponent_maps_);

                std::vector<fea_io::ArrayInfo> array_infos;
                for (auto& output_algorithm : output_algorithms)
                {
                    array_infos.push_back(output_algorithm->get_fea_io_array_info());
                }

                // Create MeshInfo from the mesh
                auto mesh_info = enriched_elements_mesh.get_fea_io_mesh_info();
                auto partition_info = enriched_elements_mesh.get_fea_io_partition_info();

                // Find number of quadrature points
                flex_mesh::GlobalIndex n_quad_points = 0;
                for (auto& e_ptr : mesh_->elements())
                {
                    if (e_ptr->has_formulation_type<finite_element::XFEMFormulation>())
                    {
                        n_quad_points += e_ptr->get_n_quad_points();
                    }
                }

                // Initialize the data writer for this output step
                writer.initialize(enriched_elements_path_stem.string(),
                    mesh_info, partition_info, n_quad_points, array_infos);

                auto multi_field_calc = output::get_element_multi_field_calc(output_algorithms);

                for (auto& output_algorithm : output_algorithms)
                {
                    output_algorithm->pre_loop_step(multi_field_calc);
                }

                // TODO move this to assembler?
                std::unique_ptr<finite_element::ElementWorkspace> element_workspace_ptr
                    = std::make_unique<finite_element::ElementWorkspace>();
                simple_math::GeneralVector elem_solution_vector;
                auto curr_enriched_element_local_index = 0;
                for (auto& e_ptr : mesh_->elements())
                {
                    if (e_ptr->has_formulation_type<finite_element::XFEMFormulation>())
                    {
                        // Provide the current solution in case the element requries it.
                        elem_solution_vector = gl_displacement_->get_element_part(e_ptr->get_local_number());

                        // This should be separated out in the future, but for now provide the workspace
                        e_ptr->set_workspace(element_workspace_ptr.get());

                        // Start calculation block
                        e_ptr->start_calculations(elem_solution_vector);

                        auto& xfem_form =
                            static_cast<finite_element::XFEMFormulation&>(
                                e_ptr->get_formulation());
                        for (auto leaf_index : utility::range_from_zero(xfem_form.get_n_leaves()))
                        {
                            xfem_form.calculate(multi_field_calc, leaf_index);

							auto poly_e_ptr = enriched_elements_mesh.elements()[curr_enriched_element_local_index];

                            for (auto& output_algorithm : output_algorithms)
                            {
                                // Use the element in the polyhedral mesh for the output algorithm
                                output_algorithm->add_element_contribution(*poly_e_ptr, multi_field_calc);
                            }

							curr_enriched_element_local_index++;
                        }

                        // End calculation block
                        e_ptr->end_calculations();

                    }
                }

                for (auto& output_algorithm : output_algorithms)
                {
                    output_algorithm->post_loop_step();
                }

                writer.close();
                ascii_writer.close();
            }
        }

        output_meta_.update_meta_file();
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

            //std::string output_path = "ElementStiffness_" + std::to_string(fe.get_global_number()) + ".txt";
            //std::ofstream output_stream(output_path);
            //output_stream << element_matrix.str() << std::endl;
            //output_stream.close();

            return element_matrix;
        });

        /*
        {
            std::string output_path = "DoFMaps.txt";
            std::ofstream output_stream(output_path);
            output_stream << dof_map_->str() << std::endl;
            output_stream.close();
        }
        */

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

		std::cout << dof_map_->str() << std::endl;

        // TODO make the solver automatically update the ghost values of the solution
        sparse_solver_->solve(
            gl_stiffness_,
            gl_rhs_,
            gl_displacement_);

        //auto eigen_vec = gl_displacement_->get_as_eigen_vector();
        //std::cout << "Solution:\n" << eigen_vec << std::endl;

        // Debug
        //auto solution = gl_displacement_->get_as_eigen_vector();
        //fea_crack::write_enrichment_surfaces_vtk_tmp(
        //    *mesh_, *dof_map_, "Solution.vtk", solution);
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
                    auto new_crack_ptr = enrichment_surfaces_->add_enrichment_surface();
                    new_crack_ptr->initiate(
                        *e_ptr,
                        0,
                        initiation_result.point,
                        initiation_result.crack_plane_normal);
                    successful = true;
                }
            }

            auto xfem_formulation = dynamic_cast<finite_element::XFEMFormulation*>(&formulation);
            if (xfem_formulation)
            {
                auto initiation_results = xfem_formulation->eval_initiation_criterion();
                for (auto i : utility::indices(initiation_results))
                {
                    auto& initiation_result = initiation_results[i];
                    if (initiation_result.criterion_satisified)
                    {
                        auto new_crack_ptr = enrichment_surfaces_->add_enrichment_surface();
                        new_crack_ptr->initiate(
                            *e_ptr,
                            i,
                            initiation_result.point,
                            initiation_result.crack_plane_normal);
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
