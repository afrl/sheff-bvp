#include "fea_analysis/FractureAnalysis.h"

#include "flex_mesh/maps/AdjacencyLists.h"
#include "flex_mesh/algorithms/Search.h"
#include "flex_mesh/element/PolyhedronElement.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/IOFeaIo.h"

#include "simple_math/Utility.h"

#include "finite_element/FiniteElement.h"
#include "finite_element/field_layouts/HierarchicalFieldLayout.h"

#include "sparse_linear_algebra/LinAlgFactory.h"

#include "output/Output.h"
#include "output/OutputFeaAnalysisHelpers.h"

#include "constitutive/enrichment/initiation_criteria/InfiniteStrengthInitiationCriterion.h"
#include "constitutive/enrichment/growth_criteria/AlwaysGrowGrowthCriterion.h"
#include "constitutive/enrichment/direction_criteria/MaxPrincipalStressDirectionCriterion.h"

#include "fea_crack/dof_algorithm/CGXFEMDofAlgorithm.h"
#include "fea_crack/misc/JsonConversion.h"

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
            fea_crack::CGXFEMDofAlgorithm::get(*dof_map_));

        dof_map_->start_enrichment_phase();
    }

    void FractureAnalysis::insert_precracks(
        const std::vector<SurfaceInitiationAtPoint>& initations)
    {
        auto e_workspace = std::make_unique<finite_element::ElementWorkspace>();

        auto& factory = finite_element::Factory::get();
        for (auto& initation : initations)
        {
            //////////////////////////////////////////////////////////////////////////////////////
            // Create and initiate the surface

            auto surface_model = factory.get_enrichment_surface_model(
                "InfinitesimalStrain",
                initation.properties,
                initation.growth_criterion,
                initation.direction_criterion);

            auto& new_surface = enrichment_surfaces_->add_enrichment_surface(
                surface_model.enrichment_model,
                surface_model.surface_model_data);

            auto initial_element_ptr = flex_mesh::find_element_with_point(*mesh_, initation.initiation_point);
            if (!initial_element_ptr)
            {
                throw _UTIL_MSG_EXCEPT("No element was found that contains the point "
                    + simple_math::to_string(initation.initiation_point));
            }

            auto& direction_criteria = initation.direction_criterion;
            simple_math::VoigtVector zero_stress = simple_math::VoigtVector::Zero(6);
            auto crack_normal = direction_criteria.eval(
                new_surface.get_global_id(),
                initation.initiation_point,
                zero_stress);

            new_surface.initiate(*initial_element_ptr, 0, initation.initiation_point, crack_normal);

            enrichment_surfaces_->write_enrichments_vtk_file(
                "./",
                new_surface.get_global_id(), // step #
                0);  // substep #
        }

        for (auto& surface : *enrichment_surfaces_)
        {
            //////////////////////////////////////////////////////////////////////////////////////
            // Grow the new surface until it cannot grow any further

            auto& surface_state = surface.get_representation().get_state();
            auto& enrichments_at_front = surface.get_enrichments_at_front();
            int growth_step = 1;
            while (!enrichments_at_front.empty())
            {
                std::vector<fea_crack::EnrichmentAtFront*> enrichments_at_front_to_grow;
                for (auto& enrichment_at_front : enrichments_at_front)
                {
                    auto& element = enrichment_at_front->get_enrichment_ref().get_element();
                    auto& criterion = enrichment_at_front->get_enrichment_ref().get().get_growth_criterion();
                    if (criterion.eval(
                            surface_state.get_enrichment_surface_id(),
                            element.get_mean_coordinate(),
                            simple_math::Vector3D::Zero()))
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
                    auto search_result = surface.search_ahead_of_front(
                        *enrichment_at_front_to_grow);

                    surface.grow_crack_front_into_all_leaves(
                        *enrichment_at_front_to_grow,
                        search_result,
                        *e_workspace);
                }
                
                enrichment_surfaces_->write_enrichments_vtk_file(
                    "./",
                    surface.get_global_id(), // step #
                    static_cast<int>(growth_step));  // substep #

                growth_step++;
            }
        }
    }

    void FractureAnalysis::insert_precracks_by_function(
        std::vector<SurfaceInsertionByLevelSet>& insertions)
    {
        auto& factory = finite_element::Factory::get();
        for (auto& insertion : insertions)
        {
            auto surface_model = factory.get_enrichment_surface_model(
                "InfinitesimalStrain",
                insertion.properties,
                constitutive::AlwaysGrowGrowthCriterion::get(),
                constitutive::MaxPrincipalStressDirectionCriterion::get());

            auto& new_surface = enrichment_surfaces_->add_enrichment_surface(
                surface_model.enrichment_model,
                surface_model.surface_model_data);

            new_surface.insert_based_on_level_set(
                insertion.level_set_function,
                [](const simple_math::Point&) { return true; },
                true,
                insertion.properties_for_fields_on_neg_side);

            enrichment_surfaces_->write_enrichments_vtk_file(
                "./",
                new_surface.get_global_id(), // step #
                0);  // substep #

        }
    }

    void FractureAnalysis::insert_precracks_by_function_outside_previous(
        const std::vector<SurfaceInsertionByLevelSet>& insertions)
    {
        auto& factory = finite_element::Factory::get();
        for (auto i : utility::indices(insertions))
        {
            auto& insertion = insertions[i];

            auto surface_model = factory.get_enrichment_surface_model(
                "InfinitesimalStrain",
                insertion.properties,
                constitutive::AlwaysGrowGrowthCriterion::get(),
                constitutive::MaxPrincipalStressDirectionCriterion::get());

            auto& new_surface = enrichment_surfaces_->add_enrichment_surface(
                surface_model.enrichment_model,
                surface_model.surface_model_data);

            new_surface.insert_based_on_level_set(
                insertion.level_set_function,
                [&insertions, &i](const simple_math::Point& p)
                {
                    for (auto j : utility::range_from_zero(i))
                    {
                        if (insertions[j].level_set_function(p) <= 0.0)
                        {
                            return false;
                        }
                    }
                    return true;
                },
                false,
                insertion.properties_for_fields_on_neg_side);

            enrichment_surfaces_->write_enrichments_vtk_file(
                "./",
                new_surface.get_global_id(), // step #
                0);  // substep #

        }
    }

    void FractureAnalysis::insert_precracks_by_function_with_boolean(
        const std::vector<SurfaceInsertionByLevelSet>& insertions)
    {
        // TODO relax this with further development.
        if (insertions.size() != 2)
        {
            // Throw error
        }

        // TODO
        auto& f1 = insertions[0].level_set_function;
        auto& f2 = insertions[1].level_set_function;
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

        auto& factory = finite_element::Factory::get();

        auto surface_model = factory.get_enrichment_surface_model(
            "InfinitesimalStrain",
            insertions[0].properties,
            constitutive::AlwaysGrowGrowthCriterion::get(),
            constitutive::MaxPrincipalStressDirectionCriterion::get());

        // Insert surface for f1 (including boolean)
        {
            auto& new_surface = enrichment_surfaces_->add_enrichment_surface(
                surface_model.enrichment_model,
                surface_model.surface_model_data);
            new_surface.insert_based_on_level_set(
                f1_wrapped,
                [](const simple_math::Point& p) { return true; });
        }

        enrichment_surfaces_->write_enrichments_vtk_file(
            "./",
            0, // step #
            0);  // substep #

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

            auto& new_surface = enrichment_surfaces_->add_enrichment_surface(
                surface_model.enrichment_model,
                surface_model.surface_model_data);
            new_surface.insert_based_on_level_set(
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

        enrichment_surfaces_->write_enrichments_vtk_file(
            "./",
            0, // step #
            1);  // substep #
    }

    void FractureAnalysis::pre_solve_tasks()
    {
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

	void FractureAnalysis::try_to_grow_precracks()
	{
        auto e_workspace = std::make_unique<finite_element::ElementWorkspace>();

		// Grow all the cracks one element until they can no longer grow
		int growth_step = 1;
		for (auto& surface : *enrichment_surfaces_)
		{
			auto& surface_state = surface.get_representation().get_state();
			auto& enrichments_at_front = surface.get_enrichments_at_front();
            auto& growth_criterion = surface.get_model().pointwise_model_data.growth_criterion;
			while (!enrichments_at_front.empty())
			{
				std::vector<fea_crack::EnrichmentAtFront*> enrichments_at_front_to_grow;
				for (auto& enrichment_at_front : enrichments_at_front)
				{
					auto& element = enrichment_at_front->get_enrichment_ref().get_element();
					if (growth_criterion.get().eval(surface_state.get_enrichment_surface_id(),
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
					auto search_result = surface.search_ahead_of_front(
						*enrichment_at_front_to_grow);

                    surface.grow_crack_front_into_all_leaves(
						*enrichment_at_front_to_grow,
						search_result,
                        *e_workspace);
				}
                surface.write_vtk("./", 1, growth_step);
				growth_step++;
			}
		}

		for (auto& crack : *enrichment_surfaces_)
		{
			crack.write_vtk("./", 1, 0);
		}

        fea_crack::write_json_enrichment_info(
            "./enrichments.json",
            mesh_->elements(),
            *dof_map_);
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
            if (e_ptr->get_n_enrichments() > 0)
            {
                finite_element::EnrichmentTreeTraverser tree_traverser(*e_ptr);

                auto& hhe_layout = dynamic_cast<const finite_element::HierarchicalFieldLayout&>(
                    e_ptr->get_field_layout());
                auto& domains = hhe_layout.get_domains();

                flex_mesh::NodeX3Matrix node_coords;
                e_ptr->get_node_coordinates(node_coords);

                // Add nodes and elements for the integration domain of each field
                // Fills enriched_elements_mesh
                for (auto leaf_index : utility::range_from_zero(e_ptr->get_n_basal_fields()))
                {
                    auto leaf_info = tree_traverser.get_leaf(leaf_index);
                    auto& field_domain = domains.get_field_domain(leaf_info.field_index);
                    auto field_polyhedron = field_domain.get_polyhedron();

                    // Add points for polyhedra
                    flex_mesh::Connectivity<Node, Element::ConnSize> nodes_for_points;
                    std::map<Element::ConnSize, Element::ConnSize> poly_to_facet_point_index_map;
                    for (auto point_index : utility::range_from_zero(field_polyhedron.points.size()))
                    {
                        // The points of the integration domain is in the parametric coordinate
                        // system, so they need to be transformed to the global for plotting.
                        auto& p = field_polyhedron.points[point_index];
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
                    }

                    std::vector<utility::DynamicArray<Element::ConnSize, Element::ConnSize>> indices_for_facets;
                    for (auto facet_index : utility::indices(field_polyhedron.point_indices_by_facet))
                    {
                        auto& facet_point_indices = field_polyhedron.point_indices_by_facet[facet_index];
                        indices_for_facets.emplace_back();
                        indices_for_facets.back().resize(facet_point_indices.size());
                        for (auto i : utility::range_from_zero(facet_point_indices.size()))
                        {
                            indices_for_facets.back()[i] = facet_point_indices[i];
                        }
                    }

                    auto poly_element = std::make_unique<flex_mesh::PolyhedronElement<Node>>();
                    poly_element->initialize(nodes_for_points, indices_for_facets);
                    poly_element->set_local_number(curr_local_element_i++);
                    poly_element->set_global_number(e_ptr->get_global_number());
                    poly_element->set_component_id(e_ptr->get_component_id());
                    enriched_elements_mesh.elements().add(std::move(poly_element));
                }

                // Add nodes and elements for each enrichment
                // Fills enrichment_surfaces_mesh
                for (auto enrichment_index : utility::range_from_zero(e_ptr->get_n_enrichments()))
                {
                    auto patches = domains.get_enrichment_domain(enrichment_index).get_patches();
                    for (auto& patch : patches)
                    {
                        std::vector<Node*> connectivity;
                        auto& points = patch.points;
                        for (auto i : utility::range_from_zero(points.size()))
                        {
                            // The points of the integration domain is in the parametric coordinate
                            // system, so they need to be transformed to the global for plotting.
                            auto& p = points[i];
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
                            flex_mesh::LCSAssignmentType::None);

                        polygon_element->set_local_number(curr_surf_local_element_i++);
                        polygon_element->set_global_number(
                            tree_traverser.get_enrichment_ref(enrichment_index).get().surface_id);
                        polygon_element->set_component_id(e_ptr->get_component_id());
                        // Right now, this throws off the reader since the component selection
                        // does not apply to enrichment surfaces.  How do we want to indicate surface?
                        //tree_traverser.get_enrichment(enrichment_index).enrichment_surface_id);
                        enrichment_surfaces_mesh.elements().add(std::move(polygon_element));
                    }
                }
            }
        }
    }

    void FractureAnalysis::write_output()
    {
        LinearElasticAnalysis::write_output();

        surface_output_flags_.push_back(
            std::make_unique<output::OutputFlag<finite_element::SurfaceQuantity>>(
                "DisplacementJump", "Global", "Node", output_flag_cache_));
        surface_output_flags_.push_back(
            std::make_unique<output::OutputFlag<finite_element::SurfaceQuantity>>(
                "EnrichmentSurfaceID", "", "Element", output_flag_cache_));

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
        enrichment_surfaces_mesh_object.data_paths.push_back(enrichment_surfaces_rel_stem);

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
                fea_io::DataWriter<flex_mesh::GlobalIndex, flex_mesh::LocalIndex, flex_mesh::FloatType> volume_writer;
                fea_io::DataWriter<flex_mesh::GlobalIndex, flex_mesh::LocalIndex, flex_mesh::FloatType> surface_writer;

                std::vector<output::OutputFlag<finite_element::VolumeQuantity>*> volume_flags_to_be_output;
                for (auto& flag : volume_output_flags_)
                {
                    if (flag->should_be_output())
                    {
                        volume_flags_to_be_output.push_back(flag.get());
                    }
                }
                std::vector<output::OutputFlag<finite_element::SurfaceQuantity>*> surface_flags_to_be_output;
                for (auto& flag : surface_output_flags_)
                {
                    if (flag->should_be_output())
                    {
                        surface_flags_to_be_output.push_back(flag.get());
                    }
                }

                // TODO only update this if the mesh has changed in subsequent writes
                flex_mesh::SubcomponentMaps enriched_elements_subcomponent_maps;
                enriched_elements_subcomponent_maps.update_maps(enriched_elements_mesh);
                flex_mesh::SubcomponentMaps enrichment_surfaces_subcomponent_maps;
                enrichment_surfaces_subcomponent_maps.update_maps(enrichment_surfaces_mesh);

                std::ofstream volume_ascii_writer((enriched_elements_path_stem / ".txt").string());
                std::ofstream surface_ascii_writer((enriched_elements_path_stem / ".txt").string());

                // Create an output algorithm for each output flag involved in this write
                auto volume_output_algorithms = output::get_output_algorithms_for_volumes(
                    volume_ascii_writer, volume_writer,
                    enriched_elements_mesh.get_partition_id(), volume_flags_to_be_output,
                    static_cast<flex_mesh::LocalIndex>(enriched_elements_mesh.nodes().size()),
                    subcomponent_maps_);
                auto surface_output_algorithms = output::get_output_algorithms_for_surfaces(
                    surface_ascii_writer, surface_writer,
                    enrichment_surfaces_mesh.get_partition_id(), surface_flags_to_be_output,
                    static_cast<flex_mesh::LocalIndex>(enrichment_surfaces_mesh.nodes().size()),
                    subcomponent_maps_);

                std::vector<fea_io::ArrayInfo> volume_array_infos;
                for (auto& output_algorithm : volume_output_algorithms)
                {
                    volume_array_infos.push_back(output_algorithm->get_fea_io_array_info());
                }
                std::vector<fea_io::ArrayInfo> surface_array_infos;
                for (auto& output_algorithm : surface_output_algorithms)
                {
                    surface_array_infos.push_back(output_algorithm->get_fea_io_array_info());
                }

                // Create MeshInfo from the mesh
                auto volume_mesh_info = enriched_elements_mesh.get_fea_io_mesh_info();
                auto volume_partition_info = enriched_elements_mesh.get_fea_io_partition_info();
                auto surface_mesh_info = enrichment_surfaces_mesh.get_fea_io_mesh_info();
                auto surface_partition_info = enrichment_surfaces_mesh.get_fea_io_partition_info();

                // Find number of quadrature points
                flex_mesh::GlobalIndex volume_n_quad_points = 0;
                for (auto& e_ptr : mesh_->elements())
                {
                    if (e_ptr->get_n_enrichments() > 0)
                    {
                        for (auto& basal_field : e_ptr->basal_fields())
                        {
                            volume_n_quad_points += basal_field.model_data.isv.size();
                        }
                    }
                }

                // Initialize the data writer for this output step
                volume_writer.initialize(enriched_elements_path_stem.string(),
                    volume_mesh_info, volume_partition_info,
                    volume_n_quad_points, volume_array_infos);
                surface_writer.initialize(enrichment_surfaces_path_stem.string(),
                    surface_mesh_info, surface_partition_info,
                    0, surface_array_infos);

                auto volume_calc_flags =
                    output::get_calc_flags_for_volumes(volume_output_algorithms);
                auto surface_calc_flags =
                    output::get_calc_flags_for_surfaces(surface_output_algorithms);
                for (auto& output_algorithm : volume_output_algorithms)
                {
                    output_algorithm->pre_loop_step(volume_calc_flags);
                }
                for (auto& output_algorithm : surface_output_algorithms)
                {
                    output_algorithm->pre_loop_step(surface_calc_flags);
                }

                // TODO move this to assembler?
                std::unique_ptr<finite_element::ElementWorkspace> element_workspace_ptr
                    = std::make_unique<finite_element::ElementWorkspace>();
                simple_math::GeneralVector elem_solution_vector;
                auto curr_enriched_element_local_index = 0;
                auto curr_enrichment_surface_local_index = 0;
                for (auto& e_ptr : mesh_->elements())
                {
                    if (e_ptr->get_n_enrichments() > 0)
                    {
                        // Provide the current solution in case the element requries it.
                        elem_solution_vector = gl_displacement_->get_element_part(e_ptr->get_local_number());

                        // Start calculation block
                        e_ptr->start_calculations(
                            *element_workspace_ptr,
                            elem_solution_vector);

                        for (auto leaf_index : utility::range_from_zero(e_ptr->get_n_basal_fields()))
                        {
                            auto results = e_ptr->calculate(
                                volume_calc_flags,
                                leaf_index,
                                *element_workspace_ptr);

							auto poly_e_ptr = enriched_elements_mesh.elements()[curr_enriched_element_local_index];

                            for (auto& output_algorithm : volume_output_algorithms)
                            {
                                // Use the element in the polyhedral mesh for the output algorithm
                                output_algorithm->add_element_contribution(*poly_e_ptr, results);
                            }

							curr_enriched_element_local_index++;
                        }
                        for (auto enrichment_index : utility::range_from_zero(e_ptr->get_n_enrichments()))
                        {
                            finite_element::EnrichmentTreeTraverser tree_traverser(*e_ptr);

                            auto& hhe_layout = dynamic_cast<const finite_element::HierarchicalFieldLayout&>(
                                e_ptr->get_field_layout());

                            auto patch_results = e_ptr->calculate(
                                surface_calc_flags,
                                enrichment_index,
                                *element_workspace_ptr);

                            for (auto patch_index : utility::indices(patch_results))
                            {
                                auto polygon_e_ptr = enrichment_surfaces_mesh.elements()[curr_enrichment_surface_local_index];

                                for (auto& output_algorithm : surface_output_algorithms)
                                {
                                    // Use the element in the polygonal mesh for the output algorithm
                                    output_algorithm->add_element_contribution(*polygon_e_ptr, patch_results[patch_index]);
                                }

                                curr_enrichment_surface_local_index++;
                            }
                        }

                        // End calculation block
                        e_ptr->end_calculations(*element_workspace_ptr);

                    }
                }

                for (auto& output_algorithm : volume_output_algorithms)
                {
                    output_algorithm->post_loop_step();
                }
                for (auto& output_algorithm : surface_output_algorithms)
                {
                    output_algorithm->post_loop_step();
                }

                volume_writer.close();
                volume_ascii_writer.close();
                surface_writer.close();
                surface_ascii_writer.close();
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
        assembler_->queue_assembly(
            finite_element::FEMatricesAndVectors::Flags::LHSMatrix,
            gl_stiffness_.get());

        // Queue assembly of the global vector containing RHS terms pertaining
        // to eigenstrains (such as thermal strain, plastic strain, etc.)
        assembler_->queue_assembly(
            finite_element::FEMatricesAndVectors::Flags::RHSVector,
            gl_eigenstrain_force_.get());

        // Execute both assemblies as one element loop
        assembler_->execute(gl_displacement_.get());

        // Total RHS = constraint RHS + applied loads + inelastic RHS (temperature/moisture terms)
        *gl_rhs_ = *gl_rhs_ + *gl_eigenstrain_force_;

        ///////////
        // Solve //
        ///////////

		//std::cout << dof_map_->str() << std::endl;

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

            // Start calculation block
            e_ptr->start_calculations(
                *element_workspace_ptr,
                elem_solution_vector);

            // TODO check initiation criterion in each field here to see if it is met
            // and if so, insert a crack.
            /* Old code as an (albiet poor) example
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
            */

            // End calculation block
            e_ptr->end_calculations(*element_workspace_ptr);
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
