#pragma once

#include "bvp/bvp/BoundaryValueProblem.h"
#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/SparseSolver.h"
#include "sparse_linear_algebra/NonLinearSolver.h"
#include "sparse_linear_algebra/Vector.h"
#include "bvp/TimeDiscretization.h"

#include <memory>
#include <vector>

namespace bvp
{

struct MateriallyNonlinearLinearAlgebra final : public BVPLinearAlgebra
{
    std::unique_ptr<sparse_linear_algebra::Vector> solution;
    std::unique_ptr<sparse_linear_algebra::Vector> solution_prev;
    std::unique_ptr<sparse_linear_algebra::Vector> eigenstrain_force;
    std::unique_ptr<sparse_linear_algebra::Vector> eigenstrain_force_prev;
    std::unique_ptr<sparse_linear_algebra::Vector> rhs;
    std::unique_ptr<sparse_linear_algebra::Vector> rhs_prev;
    std::unique_ptr<sparse_linear_algebra::SparseMatrix> stiffness;
    std::unique_ptr<sparse_linear_algebra::ConstraintMatrix> constraint_matrix;
    std::unique_ptr<sparse_linear_algebra::SparseSolver> nonlinear_solve;

    virtual std::vector<sparse_linear_algebra::Vector*> vectors() const override;
    virtual std::vector<sparse_linear_algebra::SparseMatrix*> matrices() const override;
    virtual std::vector<sparse_linear_algebra::ConstraintMatrix*> constraint_matrices() const override;

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
    

    /// @brief Solves the next time step for the boundary value problem, if there any remaining.
    virtual SolveTimeStepResult solve_next_time_step(BVPData& bvp_data) const override;

    /// @brief Returns a boolean indicating if the current state is the end of time stepping.
    ///
    /// If true is returned, solve_next_time_step() will do nothing, while if false is returned,
    /// then solve_next_time_step() will solve the next time step as expected.
    virtual bool is_end_for_time_steps(BVPData& bvp_data) const override;

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
    virtual void apply_boundary_conditions(BVPData& bvp_data, float loadCoefficient) const override;
};

struct LoadStepManagerOpts
{
  /// The maximum number of iterations allowed for any iteration loop.
  int max_iterations;
  /// step reduction if convergence not achieved
  int reduction_coeff;
  /// step increase when convergence is achieved in a number of steps
  float increase_coeff;
  /// number of converged steps to start increasing the step
  int numberConvergedStepsForIncrease;
  /// load sequence from the input file????
  std::unique_ptr<sparse_linear_algebra::Vector> load_sequence_init;
};

/// A load step manager class
class NonlinearStaticLoadStepManager 
{
public:
  /// Constructor
  /// LSopts_ - from input ???
  ///
  NonlinearStaticLoadStepManager(LoadStepManagerOpts LSopt_): LSopts_(LSopts_);
  ~NonlinearStaticLoadStepManager();
  ///                      
  int number_converged_steps=0;
  int numberConsecutiveConvergedSteps=0;
  int step_number=0; 
  float loadCoefficient=0.0;
  int loadSequenceIndex=0; 
  bool previous_step_converged = false;
  std::unique_ptr<sparse_linear_algebra::Vector> load_sequence_init_current;                                     
  /// returns next load step - returns load sequence index
  virtual current_load_step() override;
  /// loadSequenceIndex+=1; loadCoefficient inclreased by the corresponding step
  virtual increase_loadSequenceIndex();
  /// returns true if it is not over the last step
  virtual bool not_last_step_done() override;
  /// returns max_iterations  
  int get_max_iterations() override;
  /// step reduction if convergence not achieved
  float get_reduction_coeff() override;
  /// returns the step increase coefficient 
  float get_increase_coeff();
  /// returns the number of converged steps to start increasing the step
  int get_numberConvergedStepsForIncrease();
  /// returns state of the previous step - converged or not and changes the state
  bool get_previous_step_state();


protected:
  LoadStepManagerOpts LSopts_;
};

}
