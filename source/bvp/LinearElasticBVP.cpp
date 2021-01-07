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
std::vector<sparse_linear_algebra::Vector*> LinearElasticLinearAlgebra::vectors() const
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

std::vector<sparse_linear_algebra::SparseMatrix*> LinearElasticLinearAlgebra::matrices() const
{
    std::vector<sparse_linear_algebra::SparseMatrix*> non_null_matrix_ptrs;
    if (stiffness != nullptr)
    {
        non_null_matrix_ptrs.push_back(stiffness.get());
    }
    return non_null_matrix_ptrs;
}

std::vector<sparse_linear_algebra::ConstraintMatrix*>
LinearElasticLinearAlgebra::constraint_matrices() const
{
    std::vector<sparse_linear_algebra::ConstraintMatrix*> non_null_contraint_matrix_ptrs;
    if (constraint_matrix != nullptr)
    {
        non_null_contraint_matrix_ptrs.push_back(constraint_matrix.get());
    }
    return non_null_contraint_matrix_ptrs;
}

void LinearElasticBVP::init_lin_algebra(BVPData& bvp_data) const
{
    auto dof_refs = create_dof_refs(bvp_data.dof_manager());

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

} // end namespace bvp
