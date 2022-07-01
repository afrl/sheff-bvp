#include "bvp/bvp/MateriallyNonlinearBVP.h"

#include "finite_element/FiniteElement.h"
#include "sparse_linear_algebra/LinAlgFactory.h"
#include "sparse_linear_algebra/NonLinearSolver.h"
#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/Vector.h"
#include "utility/Profile.h"
#include "utility/serialization/Serialization.h"

#include <chrono>

namespace bvp
{

std::vector<sparse_linear_algebra::Vector*> MateriallyNonlinearLinearAlgebra::vectors() const
{
    std::vector<sparse_linear_algebra::Vector*> non_null_vector_ptrs;
    if (eigenstrain_force != nullptr)
    {
        non_null_vector_ptrs.push_back(eigenstrain_force.get());
    }
    if (rhs != nullptr)
    {
        non_null_vector_ptrs.push_back(rhs.get());
    }
    return non_null_vector_ptrs;
}

std::vector<sparse_linear_algebra::SparseMatrix*> MateriallyNonlinearLinearAlgebra::matrices() const
{
    std::vector<sparse_linear_algebra::SparseMatrix*> non_null_matrix_ptrs;
    if (stiffness != nullptr)
    {
        non_null_matrix_ptrs.push_back(stiffness.get());
    }
    return non_null_matrix_ptrs;
}

std::vector<sparse_linear_algebra::ConstraintMatrix*>
MateriallyNonlinearLinearAlgebra::constraint_matrices() const
{
    std::vector<sparse_linear_algebra::ConstraintMatrix*> non_null_contraint_matrix_ptrs;
    if (constraint_matrix != nullptr)
    {
        non_null_contraint_matrix_ptrs.push_back(constraint_matrix.get());
    }
    return non_null_contraint_matrix_ptrs;
}

void MateriallyNonlinearBVP::init_lin_algebra(BVPData& bvp_data) const
{
    _F_BEGIN2("allocate global vectors/matrices");

    auto dof_refs = create_dof_refs(bvp_data.dof_manager());

    bvp_data.gl_lin_alg = std::make_unique<MateriallyNonlinearLinearAlgebra>();
    auto& lin_alg       = static_cast<MateriallyNonlinearLinearAlgebra&>(*bvp_data.gl_lin_alg);

    // Right now, penalty is all that works.  Switch this to in-place elimination when that works.
    lin_alg.constraint_matrix =
      sparse_linear_algebra::LinAlgFactory::get().create_constraint_matrix(
        dof_refs,
        bvp_data.core_data.solver_opts.constraint_method);

    lin_alg.stiffness = sparse_linear_algebra::LinAlgFactory::get().create_matrix(dof_refs);
    lin_alg.stiffness->set_constraint_matrix(lin_alg.constraint_matrix.get());
    lin_alg.stiffness->allocate();

    lin_alg.delta_solution = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.delta_solution->allocate();
    lin_alg.solution = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.solution->allocate();
    lin_alg.solution_prev = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.solution_prev->allocate();
    lin_alg.eigenstrain_force = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.eigenstrain_force->allocate();
    lin_alg.eigenstrain_force_prev =
      sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.eigenstrain_force_prev->allocate();
    lin_alg.rhs = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.rhs->allocate();
    lin_alg.rhs_prev = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.rhs_prev->allocate();

    lin_alg.sparse_solver =
      sparse_linear_algebra::LinAlgFactory::get().create_solver(bvp_data.core_data.solver_opts);

    _F_END();
}

void MateriallyNonlinearBVP::apply_boundary_conditions(BVPData& bvp_data) const
{
    auto& lin_alg = static_cast<MateriallyNonlinearLinearAlgebra&>(*bvp_data.gl_lin_alg);

    for (auto& bc : bvp_data.bcs())
    {
        const auto scale_factor_function_index = bc->get_scale_factor_function_index();
        auto bc_scale_factor                   = 0.0;
        if (scale_factor_function_index != -1)
        {
            const auto current_time = bvp_data.core_data.time_discr_data.current_time;
            bc_scale_factor =
              bvp_data.core_data.time_discr_data.scale_functions[scale_factor_function_index]
                .get_scale_factor(current_time);
        }

        // IMPORTANT NOTE: this assumes the field has # of components as dimensionality
        // This works for displacement with solid elements, but reconsider what should
        // be done to handle scalar fields, shell elements, combined fields, etc.
        // TODO figure out a robust term strategy for this
        bc->initialize(bvp_data.fe_mesh(),
                       bvp_data.node_to_element_map(),
                       bvp_data.dof_manager().get_dof_maps(),
                       utility::GlobalVariables::get_dimensionality(),
                       bc_scale_factor,
                       true, /* is_incremental */
                       *lin_alg.constraint_matrix);
    }
}

SolveTimeStepResult MateriallyNonlinearBVP::solve_next_time_step(BVPData& bvp_data) const
{
    init_all_data(bvp_data);

    auto& time_discr_data = bvp_data.core_data.time_discr_data;

    auto current_time_step   = time_discr_data.last_converged_time_step + 1;
    double target_time       = time_discr_data.target_times[current_time_step];
    double last_substep_time = time_discr_data.last_converged_time;

    // Create a directory for the current time step
    auto ts_dir = bvp_data.output_dir() / ("ts_" + std::to_string(current_time_step));
    if (!fs::exists(ts_dir))
    {
        fs::create_directory(ts_dir);
    }

    auto log_stage_name = bvp_data.bvp_name() + "_ts" + std::to_string(current_time_step);
    sparse_linear_algebra::LinAlgFactory::get().start_log_stage(log_stage_name);

    auto& lin_alg = static_cast<MateriallyNonlinearLinearAlgebra&>(*bvp_data.gl_lin_alg);

    // TODO reevaluate how often the local array is updated and see if the tangent and residual can
    //      calculated in the same assembler execute call.
    auto& assembler = bvp_data.assembler();
    auto& dof_maps  = bvp_data.dof_manager().get_dof_maps();
    std::vector<double> local_solution(dof_maps.n_dofs_on_rank(), 0.0);
    std::function<void(const sparse_linear_algebra::Vector&, sparse_linear_algebra::SparseMatrix&)>
      tangent_matrix_functor =
        [&local_solution, &assembler, &dof_maps](const sparse_linear_algebra::Vector& solution,
                                                 sparse_linear_algebra::SparseMatrix& A)
    {
        solution.update_local_array(local_solution);
        assembler.queue_assembly(finite_element::FEMatricesAndVectors::Flags::LHSMatrix, &A);
        assembler.execute(local_solution, dof_maps);
    };
    std::function<void(const sparse_linear_algebra::Vector&, sparse_linear_algebra::Vector&)>
      residual_functor =
        [&local_solution, &assembler, &dof_maps](const sparse_linear_algebra::Vector& solution,
                                                 sparse_linear_algebra::Vector& residual)
    {
        solution.update_local_array(local_solution);
        // Note: this calculates the "internal force vector"
        assembler.queue_assembly(finite_element::FEMatricesAndVectors::Flags::CompleteRHSVector,
                                 &residual);
        assembler.execute(local_solution, dof_maps);
        // TODO implement residual = "external force vector" - "internal force vector"
    };

    // Note that these are references
    double& delta_time   = time_discr_data.current_delta_time;
    double& current_time = time_discr_data.current_time;
    while (last_substep_time == target_time)
    {
        // Figure out the time that will be used for this substep
        current_time = last_substep_time + delta_time;
        if (current_time > target_time)
        {
            current_time = target_time;
        }

        // Initialize boundary conditions
        apply_boundary_conditions(bvp_data);

        // try to solve nonlinear problem
        bool substep_converged = false;
        sparse_linear_algebra::NonLinearSolverParameters nl_params; // TODO set this via input file
        auto status       = sparse_linear_algebra::nonlinear_solve(nl_params,
                                                             lin_alg.solution_prev.get(),
                                                             *lin_alg.sparse_solver,
                                                             lin_alg.stiffness.get(),
                                                             lin_alg.rhs.get(),
                                                             lin_alg.delta_solution.get(),
                                                             tangent_matrix_functor,
                                                             residual_functor);
        substep_converged = status.integral_stats["converged"] == 1;

        if (substep_converged)
        {
            last_substep_time = current_time;

            *lin_alg.solution_prev = *lin_alg.solution;
            *lin_alg.solution      = *lin_alg.delta_solution + *lin_alg.solution_prev;

            /// update the internal state variables (overwritting the
            /// internal state variables for the last converged time step)
            assembler.queue_operation(
              [](finite_element::FiniteElement::ConstParams& params,
                 finite_element::FiniteElementAccessor& fe)
              { finite_element::FiniteElement::update_converged_isvs(params, fe.layout); });
            lin_alg.solution->update_local_array(local_solution);
            // Execute the update for ISVs
            assembler.execute(local_solution, dof_maps);

            if (time_discr_data.n_consecutive_converged_substeps ==
                time_discr_data.n_consecutive_substeps_before_increase)
            {
                time_discr_data.n_consecutive_converged_substeps = 0;
                delta_time *= 2.0;
            }
        }
        else
        {
            time_discr_data.n_consecutive_converged_substeps = 0;
            delta_time /= 2.0;
        }
    }

    sparse_linear_algebra::LinAlgFactory::get().end_current_log_stage();

    return { true, /* converged */
             time_discr_data.current_time ==
               time_discr_data.target_times.back() /* done with time steps */ };
}

} // end namespace bvp
