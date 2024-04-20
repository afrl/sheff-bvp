#pragma once

#include "bvp/bvp/BoundaryValueProblem.h"
#include "bvp/bvp/FractureBVP.h"
#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/SparseSolver.h"
#include "sparse_linear_algebra/Vector.h"

#include <memory>
#include <vector>

namespace bvp
{

struct MateriallyNonlinearLinearAlgebra final : public BVPLinearAlgebra
{
    std::unique_ptr<sparse_linear_algebra::Vector> delta_solution;
    std::unique_ptr<sparse_linear_algebra::Vector> solution;
    std::unique_ptr<sparse_linear_algebra::Vector> solution_prev;
    std::unique_ptr<sparse_linear_algebra::Vector> eigenstrain_force;
    std::unique_ptr<sparse_linear_algebra::Vector> eigenstrain_force_prev;
    std::unique_ptr<sparse_linear_algebra::Vector> rhs;
    std::unique_ptr<sparse_linear_algebra::Vector> rhs_prev;
    std::unique_ptr<sparse_linear_algebra::Vector> residual;
    std::unique_ptr<sparse_linear_algebra::SparseMatrix> stiffness;
    std::unique_ptr<sparse_linear_algebra::ConstraintMatrix> constraint_matrix;
    std::unique_ptr<sparse_linear_algebra::SparseSolver> sparse_solver;

    virtual std::vector<sparse_linear_algebra::Vector*> vectors() const override;
    virtual std::vector<sparse_linear_algebra::SparseMatrix*> matrices() const override;
    virtual std::vector<sparse_linear_algebra::ConstraintMatrix*> constraint_matrices()
      const override;
};

/// A materially nonlinear analysis
class MateriallyNonlinearBVP : public BoundaryValueProblem
{
  public:
    static MateriallyNonlinearBVP& get()
    {
        static MateriallyNonlinearBVP singleton;
        return singleton;
    }

    virtual ~MateriallyNonlinearBVP() {}

    /// TODO fix this...does FractureBVP need to inherit from this?
    /// Insert precracks by initiating a new crack at the given point and growing the surface
    /// via the direction criterion and growth criterion that were assigned to the elements.
    void insert_precracks(BVPData& bvp_data,
                          const std::vector<SurfaceInitiationAtPoint>& initations) const;

    /// Insert precracks along implicit functions that define the surface (where the functions
    /// equals 0).  This method differs from insert_precracks_by_function in that it only
    /// inserts a surface where the sign of previously inserted surfaces (in the vector) are
    /// positive.  For example, if the first function defines a sphere with the normal pointing
    /// outward and the second function defines a plane intersecting the sphere, then a surface
    /// will only be inserted outside the sphere.
    void insert_precracks_by_function_outside_previous(
      BVPData& bvp_data,
      const std::vector<SurfaceInsertionByLevelSet>& insertions) const;


    /// @brief Solves the next time step for the boundary value problem, if there any remaining.
    virtual SolveTimeStepResult solve_next_time_step(BVPData& bvp_data) const override;

  protected:
    MateriallyNonlinearBVP() {}

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

    /// @brief If enrichment has occurred, this renumbers the DoFs to satisfy the needs of the
    /// solvers and enforce continuity across rank boundaries.
    void renumber_dofs(BVPData& bvp_data) const;

    // TODO - again deal with combining FractureBVP
    xfem::Params get_xfem_params(BVPData& bvp_data) const;
};

}
