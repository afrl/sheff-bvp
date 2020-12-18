#include "bvp/bvp/LinearElasticBVP.h"

#include "finite_element/FiniteElement.h"
#include "sparse_linear_algebra/LinAlgFactory.h"
#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/Vector.h"
#include "utility/Profile.h"
#include "utility/serialization/Serialization.h"

#include <chrono>

namespace bvp
{
void LinearElasticBVP::init_lin_algebra(BVPData& bvp_data) const
{
    auto& dof_maps                             = bvp_data.dof_manager().get_dof_maps();
    auto& enrichment_data                      = bvp_data.dof_manager().get_enrichment_data();
    sparse_linear_algebra::DofMapRefs dof_refs = {
        dof_maps.n_owned_dofs,           dof_maps.n_local_ghost_dofs,
        dof_maps.n_exclusive_ghost_dofs, dof_maps.owned_global_dof_begin,
        dof_maps.owned_global_dof_end,   dof_maps.n_owned_elements,
        dof_maps.element_to_rank_maps,   dof_maps.element_n_dofs,
        dof_maps.rank_to_global_map,     enrichment_data.new_dof_to_old_rank_map
    };

    bvp_data.gl_lin_alg = std::make_unique<LinearElasticLinearAlgebra>();
    auto& lin_alg       = static_cast<LinearElasticLinearAlgebra&>(*bvp_data.gl_lin_alg);

    // Right now, penalty is all that works.  Switch this to in-place elimination when that works.
    lin_alg.constraint_matrix = sparse_linear_algebra::LinAlgFactory::create_constraint_matrix(
      dof_refs,
      sparse_linear_algebra::ConstraintMethod::Penalty);

    lin_alg.stiffness = sparse_linear_algebra::LinAlgFactory::create_matrix(dof_refs);
    lin_alg.stiffness->set_constraint_matrix(lin_alg.constraint_matrix.get());
    lin_alg.stiffness->allocate();

    lin_alg.eigenstrain_force = sparse_linear_algebra::LinAlgFactory::create_vector(dof_refs);
    lin_alg.eigenstrain_force->allocate();
    lin_alg.rhs = sparse_linear_algebra::LinAlgFactory::create_vector(dof_refs);
    lin_alg.rhs->allocate();

    sparse_linear_algebra::SparseSolverOptions opts;
    opts.relative_tol          = 1e-16;
    opts.requires_nonsymmetric = false;
    lin_alg.sparse_solver      = sparse_linear_algebra::LinAlgFactory::create_solver(opts);
}

void LinearElasticBVP::apply_boundary_conditions(BVPData& bvp_data) const
{
    auto& lin_alg = static_cast<LinearElasticLinearAlgebra&>(*bvp_data.gl_lin_alg);

    for (auto& bc : bvp_data.bcs())
    {
        // IMPORTANT NOTE: this assumes the field has # of components as dimensionality
        // This works for displacement with solid elements, but reconsider what should
        // be done to handle scalar fields, shell elements, combined fields, etc.
        // TODO figure out a robust term strategy for this
        bc->initialize(bvp_data.fe_mesh(),
                       bvp_data.node_to_element_map(),
                       bvp_data.dof_manager().get_dof_maps(),
                       utility::GlobalVariables::get_dimensionality(),
                       *lin_alg.constraint_matrix);
    }
}

SolveTimeStepResult LinearElasticBVP::solve_next_time_step(BVPData& bvp_data) const
{
    _F_BEGIN();

    init_all_data(bvp_data);

    // Write restart data
    // Create a directory for the current time step
    bvp_data.core_data.time_discr_data.last_converged_time_step++;
    auto curr_time_step = bvp_data.core_data.time_discr_data.last_converged_time_step;
    auto ts_dir         = bvp_data.output_dir() / ("ts_" + std::to_string(curr_time_step));
    if (!fs::exists(ts_dir))
    {
        fs::create_directory(ts_dir);
    }

    _F_BEGIN2("apply BC's");

    // TODO indicate if Secant or Tangent with BCUpdateType::Secant
    apply_boundary_conditions(bvp_data);

    _F_END();

    //////////////
    // Assembly //
    //////////////

    _F_BEGIN2("assembly");

    auto& assembler            = bvp_data.assembler();
    auto& gl_displacement      = bvp_data.solution();
    auto& lin_alg              = static_cast<LinearElasticLinearAlgebra&>(*bvp_data.gl_lin_alg);
    auto& gl_stiffness         = *lin_alg.stiffness;
    auto& gl_eigenstrain_force = *lin_alg.eigenstrain_force;
    auto& gl_rhs               = *lin_alg.rhs;

    // Queue assembly of the global stiffness matrix
    assembler.queue_assembly(finite_element::FEMatricesAndVectors::Flags::LHSMatrix, &gl_stiffness);

    // Queue assembly of the global vector containing RHS terms pertaining
    // to eigenstrains (such as thermal strain, plastic strain, etc.)
    assembler.queue_assembly(finite_element::FEMatricesAndVectors::Flags::RHSVector,
                             &gl_eigenstrain_force);

    // Execute both assemblies as one element loop
    assembler.execute(&gl_displacement);

    _F_END();

    // Total RHS = constraint RHS + applied loads + inelastic RHS (temperature/moisture terms)
    gl_rhs = gl_rhs + gl_eigenstrain_force;

    if (BVPGlobalOptions::get().save_matrices)
    {
        gl_stiffness.save((ts_dir / "matrix_pre_constraints.bin").generic_string());
        gl_rhs.save((ts_dir / "rhs_pre_constraints.bin").generic_string());
    }

    ///////////
    // Solve //
    ///////////

    _F_BEGIN2("linear algebra solve");

    // TODO make the solver automatically update the ghost values of the solution
    lin_alg.sparse_solver->solve(&gl_stiffness, &gl_rhs, &gl_displacement);

    _F_END();
    
    if (BVPGlobalOptions::get().save_matrices)
    {
        gl_stiffness.save((ts_dir / "matrix_post_constraints.bin").generic_string());
        gl_rhs.save((ts_dir / "rhs_post_constraints.bin").generic_string());
        gl_displacement.save((ts_dir / "solution.bin").generic_string());
    }

    save_bvp_data(bvp_data, (ts_dir / "converged").generic_string());

    _F_END();

    return { true, is_end_for_time_steps(bvp_data) };
}

bool LinearElasticBVP::is_end_for_time_steps(BVPData& bvp_data) const
{
    auto& time_discr_data = bvp_data.core_data.time_discr_data;
    if (time_discr_data.type == TimeDiscretizationType::STATIC)
    {
        return time_discr_data.last_converged_time_step == 0;
    }
    throw _UTIL_MSG_EXCEPT_C("Discretization type not implimented.");
}

/*
*
void LinearElasticAnalysis::write_output()
{
    _F_BEGIN();

    auto& output_meta_info = output_meta_.get_analysis_meta();
    auto& time_step_info   = output_meta_.get_time_step(time_step_manager_->get_time_step(),
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
        time_step_info.levels_of_detail.push_back(
          { 0,                               // level
            fea_io::LODParentRelation::None, // relationship_with_parent
            { mesh_object } });              // mesh_objects

        // Create a writer for this time step
        fea_io::DataWriter<flex_mesh::GlobalIndex, flex_mesh::LocalIndex, flex_mesh::FloatType>
          writer;

        // DEBUG TODO delete
        volume_output_flags_.push_back(
          std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>("Position",
                                                                               "Global",
                                                                               "Node",
                                                                               output_flag_cache_));
        volume_output_flags_.push_back(
          std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>("Displacement",
                                                                               "Global",
                                                                               "Node",
                                                                               output_flag_cache_));
        volume_output_flags_.push_back(
          std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>("Stress",
                                                                               "Local",
                                                                               "Node",
                                                                               output_flag_cache_));
        volume_output_flags_.push_back(
          std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>("Strain",
                                                                               "Global",
                                                                               "Node",
                                                                               output_flag_cache_));
        volume_output_flags_.push_back(
          std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>("Stress",
                                                                               "Global",
                                                                               "Element",
                                                                               output_flag_cache_));
        volume_output_flags_.push_back(
          std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>("Strain",
                                                                               "Global",
                                                                               "Element",
                                                                               output_flag_cache_));
        volume_output_flags_.push_back(
          std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>("Position",
                                                                               "Global",
                                                                               "QuadraturePoint",
                                                                               output_flag_cache_));
        volume_output_flags_.push_back(
          std::make_unique<output::OutputFlag<finite_element::VolumeQuantity>>("Stress",
                                                                               "Global",
                                                                               "QuadraturePoint",
                                                                               output_flag_cache_));

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
        auto output_algorithms = output::get_output_algorithms_for_volumes(
          ascii_writer_,
          writer,
          mesh_->get_partition_id(),
          flags_to_be_output,
          static_cast<flex_mesh::LocalIndex>(mesh_->nodes().size()),
          subcomponent_maps_);

        std::vector<fea_io::ArrayInfo> array_infos;
        for (auto& output_algorithm : output_algorithms)
        {
            array_infos.push_back(output_algorithm->get_fea_io_array_info());
        }

        // Create MeshInfo from the mesh
        auto mesh_info      = mesh_->get_fea_io_mesh_info();
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
                          mesh_info,
                          partition_info,
                          n_quad_points,
                          array_infos);

        auto volume_calc_flags = output::get_calc_flags_for_volumes(output_algorithms);

        for (auto& output_algorithm : output_algorithms)
        {
            output_algorithm->pre_loop_step(volume_calc_flags);
        }

        // TODO move this to assembler?
        std::unique_ptr<finite_element::ElementWorkspace> element_workspace_ptr =
          std::make_unique<finite_element::ElementWorkspace>();
        simple_math::GeneralVector elem_solution_vector;
        for (auto& e_ptr : mesh_->elements())
        {
            // Provide the current solution in case the element requries it.
            elem_solution_vector = gl_displacement_->get_element_part(e_ptr->get_local_number());

            // Start calculation block
            e_ptr->start_calculations(*element_workspace_ptr, elem_solution_vector);

            auto e_number = e_ptr->get_global_number();
            auto results  = e_ptr->calculate(volume_calc_flags, *element_workspace_ptr);

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

    _F_END();
}
*/
} // end namespace bvp
