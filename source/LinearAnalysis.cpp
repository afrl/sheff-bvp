#include "fea_analysis/LinearAnalysis.h"

#include "finite_element/FiniteElement.h"

#include "sparse_linear_algebra/LinAlgFactory.h"

#include "output/Output.h"
#include "output/OutputFeaAnalysisHelpers.h"

namespace fea_analysis
{
	void LinearElasticAnalysis::initialize_linear_algebra()
    {
        dof_refs_ = std::make_unique<sparse_linear_algebra::DofMapReferences>(
            dof_map_->get_n_owned_dofs(),
            dof_map_->get_n_ghost_dofs_from_lower_ranks(),
            dof_map_->get_n_ghost_dofs_from_higher_ranks(),
            dof_map_->get_owned_global_dof_begin(),
            dof_map_->get_owned_global_dof_end(),
            dof_map_->get_n_owned_elements(),
            dof_map_->get_element_to_rank_maps(),
            dof_map_->get_element_n_dofs(),
            dof_map_->get_node_to_rank_maps(),
            dof_map_->get_node_n_dofs(),
            dof_map_->get_rank_to_global_map(),
            dof_map_->get_enrichment_new_dof_to_old_rank_map());

		// Right now, penalty is all that works.  Switch this to in-place elimination when that works.
		gl_constraint_matrix_ = sparse_linear_algebra::LinAlgFactory::create_constraint_matrix(
            *dof_refs_,
			sparse_linear_algebra::ConstraintMethod::Penalty);

        gl_stiffness_ = sparse_linear_algebra::LinAlgFactory::create_matrix(*dof_refs_);
		gl_stiffness_->set_constraint_matrix(gl_constraint_matrix_.get());
        gl_stiffness_->allocate();

        gl_displacement_ = sparse_linear_algebra::LinAlgFactory::create_vector(*dof_refs_);
        gl_displacement_->allocate();
        gl_eigenstrain_force_ = sparse_linear_algebra::LinAlgFactory::create_vector(*dof_refs_);
        gl_eigenstrain_force_->allocate();
        gl_rhs_ = sparse_linear_algebra::LinAlgFactory::create_vector(*dof_refs_);
        gl_rhs_->allocate();

        sparse_linear_algebra::SparseSolverOptions opts;
        opts.relative_tol = 1e-16;
        opts.requires_nonsymmetric = false;
        sparse_solver_ = sparse_linear_algebra::LinAlgFactory::create_solver(opts);
    }

    void LinearElasticAnalysis::solve_time_step()
    {
        // TODO indicate if Secant or Tangent with BCUpdateType::Secant
        apply_boundary_conditions();

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

        // TODO make the solver automatically update the ghost values of the solution
        sparse_solver_->solve(
            gl_stiffness_,
            gl_rhs_,
            gl_displacement_);

        /*
        auto solution = gl_displacement_->get_as_eigen_vector();
        if (utility::MPIUtils::get_local_rank() == 0)
        {
            std::ofstream out("solution.txt");
            out << solution;
            out.close();
        }
        */
    }



    void LinearElasticAnalysis::write_output()
    {
        auto& output_meta_info = output_meta_.get_analysis_meta();
        auto& time_step_info = output_meta_.get_time_step(
            time_step_manager_->get_time_step(),
            time_step_manager_->get_current_time());

        {
            fea_io::MeshObject mesh_object;
            mesh_object.type = fea_io::MeshObjectType::Volume;
            // Note: we can assign the same paths from the input info because
            // FeaAnalysis::set_input() copies all of the files over to the output directory.
            mesh_object.mesh_path =
                output_meta_info.input_info.levels_of_detail[0].mesh_objects[0].mesh_path;
            mesh_object.group_paths =
                output_meta_info.input_info.levels_of_detail[0].mesh_objects[0].group_paths;
            // Use the mesh file name for the name of the data that belongs to the mesh
            std::string mesh_filename = fs::path(mesh_object.mesh_path).filename().string();
            // Note: there is an assumption here that the time step directory is the same folder
            // as the output meta file.  This may need to be generalized.
            mesh_object.data_paths.push_back(
                (curr_time_step_dir_.filename() / mesh_filename).generic_string());
            fs::path output_data_path = curr_time_step_dir_ / mesh_filename;
            // Update the time step info
            time_step_info.levels_of_detail.push_back({
                0, // level
                fea_io::LODParentRelation::None, // relationship_with_parent
                { mesh_object } }); // mesh_objects

            // Create a writer for this time step
            fea_io::DataWriter<flex_mesh::GlobalIndex, flex_mesh::LocalIndex, flex_mesh::FloatType> writer;

            // DEBUG TODO delete
            volume_output_flags_.push_back(
                std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>(
                    "Displacement", "Global", "Node", output_flag_cache_));
            volume_output_flags_.push_back(
                std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>(
                    "Stress", "Global", "Element", output_flag_cache_));
            volume_output_flags_.push_back(
                std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>(
                    "Strain", "Global", "Element", output_flag_cache_));

            std::vector<output::OutputFlag<finite_element::VolumeQuantity>*> flags_to_be_output;
            for (auto& flag : volume_output_flags_)
            {
                if (flag->should_be_output())
                {
                    flags_to_be_output.push_back(flag.get());
                }
            }

            // TODO only update this if the mesh has changed in subsequent writes
            subcomponent_maps_.update_maps(*mesh_);

            ascii_writer_.open((output_data_path / ".txt").string());

            // Create an output algorithm for each output flag involved in this write
            auto output_algorithms = output::get_output_algorithms_for_volumes(ascii_writer_, writer,
                mesh_->get_partition_id(), flags_to_be_output,
                static_cast<flex_mesh::LocalIndex>(mesh_->nodes().size()), subcomponent_maps_);

            std::vector<fea_io::ArrayInfo> array_infos;
            for (auto& output_algorithm : output_algorithms)
            {
                array_infos.push_back(output_algorithm->get_fea_io_array_info());
            }

            // Create MeshInfo from the mesh
            auto mesh_info = mesh_->get_fea_io_mesh_info();
            auto partition_info = mesh_->get_fea_io_partition_info();

            // Find number of quadrature points
            flex_mesh::GlobalIndex n_quad_points = 0;
            for (auto& e_ptr : mesh_->elements())
            {
                for (auto& basal_field : e_ptr->basal_fields())
                {
                    n_quad_points += basal_field.model_data.isv.size();
                }
            }

            // Initialize the data writer for this output step
            writer.initialize(output_data_path.string(),
                mesh_info, partition_info, n_quad_points, array_infos);

            auto volume_calc_flags = output::get_calc_flags_for_volumes(output_algorithms);

            for (auto& output_algorithm : output_algorithms)
            {
                output_algorithm->pre_loop_step(volume_calc_flags);
            }

            // TODO move this to assembler?
            std::unique_ptr<finite_element::ElementWorkspace> element_workspace_ptr
                = std::make_unique<finite_element::ElementWorkspace>();
            simple_math::GeneralVector elem_solution_vector;
            for (auto& e_ptr : mesh_->elements())
            {
                // Provide the current solution in case the element requries it.
                elem_solution_vector = gl_displacement_->get_element_part(e_ptr->get_local_number());

                // Start calculation block
                e_ptr->start_calculations(
                    *element_workspace_ptr,
                    elem_solution_vector);

                auto e_number = e_ptr->get_global_number();
                auto results = e_ptr->calculate(volume_calc_flags, *element_workspace_ptr);

                // End calculation block
                e_ptr->end_calculations(*element_workspace_ptr);

                for (auto& output_algorithm : output_algorithms)
                {
                    output_algorithm->add_element_contribution(*e_ptr, results);
                }
            }

            for (auto& output_algorithm : output_algorithms)
            {
                output_algorithm->post_loop_step();
            }

            writer.close();
            ascii_writer_.close();
        }

        output_meta_.update_meta_file();
    }
} // end namespace fea_analysis
