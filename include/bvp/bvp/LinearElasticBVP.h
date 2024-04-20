#pragma once

#include "bvp/bvp/BoundaryValueProblem.h"
#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/SparseSolver.h"
#include "sparse_linear_algebra/Vector.h"

#include <memory>
#include <vector>

namespace bvp
{

struct LinearElasticLinearAlgebra final : public BVPLinearAlgebra
{
    std::unique_ptr<sparse_linear_algebra::Vector> solution;
    std::unique_ptr<sparse_linear_algebra::Vector> eigenstrain_force;
    std::unique_ptr<sparse_linear_algebra::Vector> rhs;
    std::unique_ptr<sparse_linear_algebra::SparseMatrix> stiffness;
    std::unique_ptr<sparse_linear_algebra::ConstraintMatrix> constraint_matrix;
    std::unique_ptr<sparse_linear_algebra::SparseSolver> sparse_solver;

    virtual std::vector<sparse_linear_algebra::Vector*> vectors() const override;
    virtual std::vector<sparse_linear_algebra::SparseMatrix*> matrices() const override;
    virtual std::vector<sparse_linear_algebra::ConstraintMatrix*> constraint_matrices()
      const override;
};

/// A linear elastic analysis
class LinearElasticBVP : public BoundaryValueProblem
{
  public:
    static LinearElasticBVP& get()
    {
        static LinearElasticBVP singleton;
        return singleton;
    }

    virtual ~LinearElasticBVP() {}


    /// @brief Solves the next time step for the boundary value problem, if there any remaining.
    virtual SolveTimeStepResult solve_next_time_step(BVPData& bvp_data) const override;

  protected:
    LinearElasticBVP() {}

    /// @brief Initializes the linear algebra containers (bvp_data.core_data.gl_solution and
    /// bvp_data.gl_lin_alg).
    ///
    /// Sets bvp_data.gl_lin_alg.initialized to true.
    virtual void init_lin_algebra(BVPData& bvp_data) const override;

    /// @brief Applies the boundary conditions in bvp_data.boundary_conditions to the system
    /// of equations
    ///
    /// Modifies the constraint matrix and appropriate vectors in bvp_data.gl_lin_alg.
    virtual void apply_boundary_conditions(BVPData& bvp_data) const override;

    /// @see BoundaryValueProblem::ensure_matrices_are_allocated
    virtual void ensure_matrices_are_allocated(BVPData& bvp_data) const override;
};
}
