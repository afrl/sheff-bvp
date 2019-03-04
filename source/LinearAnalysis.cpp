#include "fea_analysis/LinearAnalysis.h"

#include "finite_element/FiniteElement.h"
#include "finite_element/BaseFormulation.h"
#include "finite_element/elasticity/InfinitesimalStrain.h"

#include "sparse_linear_algebra/Factory.h"

namespace fea_analysis
{
	void LinearElasticAnalysis::set_element_formulations()
	{
		// For a linear analysis, for now, assume all elements get an InfinitesmalStrain formulation
		for (auto& e_ptr : mesh_->elements())
		{
			auto formulation = std::unique_ptr<finite_element::BaseFormulation, finite_element::BaseFormulationDeleter>(new finite_element::InfinitesimalStrainFormulation(e_ptr.get()));
			e_ptr->set_formulation(std::move(formulation));
			e_ptr->allocate_isvs();
		}
	}

	void LinearElasticAnalysis::initialize_linear_algebra()
    {
        sparse_linear_algebra::LinAlgFactory::init(
            sparse_linear_algebra::LinAlgFactory::LinAlgBackend::MKL);

		// Right now, penalty is all that works.  Switch this to in-place elimination when that works.
		gl_constraint_matrix_ = sparse_linear_algebra::LinAlgFactory::create_constraint_matrix(
			sparse_linear_algebra::ConstraintMethod::Penalty);

        gl_stiffness_ = sparse_linear_algebra::LinAlgFactory::create_matrix();
		gl_stiffness_->set_constraint_matrix(gl_constraint_matrix_.get());

        gl_displacement_ = sparse_linear_algebra::LinAlgFactory::create_vector();
        gl_eigenstrain_force_ = sparse_linear_algebra::LinAlgFactory::create_vector();
        gl_rhs_ = sparse_linear_algebra::LinAlgFactory::create_vector();

        initialize_container(*gl_stiffness_);
        initialize_container(*gl_displacement_);
        initialize_container(*gl_eigenstrain_force_);
        initialize_container(*gl_rhs_);

        sparse_linear_algebra::SparseSolverOptions opts;
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