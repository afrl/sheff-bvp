#include "fea_analysis/LinearAnalysis.h"

#include "finite_element/FiniteElement.h"
#include "finite_element/BaseFormulation.h"
#include "finite_element/elasticity/InfinitesimalStrain.h"

#include "sparse_linear_algebra/Factory.h"

#include "output/Output.h"
#include "output/OutputFeaAnalysisHelpers.h"

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



    void LinearElasticAnalysis::write_output()
    {
        // DEBUG TODO delete
        output::OutputFlagCache cache;
        output_flags_.push_back(std::make_unique<output::OutputFlag>("Displacement", "Global", "Node", cache));

        std::vector<output::OutputFlag*> flags_to_be_output;
        for (auto& flag : output_flags_)
        {
            if (flag->should_be_output())
            {
                flags_to_be_output.push_back(flag.get());
            }
        }

        // TODO only update this if the mesh has changed in subsequent writes
        subcomponent_maps_.update_maps(*mesh_);

        ascii_writer_.open(output::get_data_file_name() + ".txt");

        // Create an output algorithm for each output flag involved in this write
        auto output_algorithms = output::get_output_algorithms(ascii_writer_, writer_,
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
            n_quad_points += e_ptr->get_n_quad_points();
        }

        // Initialize the data writer for this output step
        writer_.initialize(output::get_data_file_name(),
            mesh_info, partition_info, n_quad_points, array_infos);

        auto multi_field_calc = output::get_element_multi_field_calc(output_algorithms);

        for (auto& output_algorithm : output_algorithms)
        {
            output_algorithm->pre_loop_step(multi_field_calc);
        }

        // TODO move this to assembler?
        finite_element::ElementWorkspace element_workspace;
        simple_math::GeneralVector elem_solution_vector;
        for (auto& e_ptr : mesh_->elements())
        {
            // Provide the current solution in case the element requries it.
            elem_solution_vector = gl_displacement_->get_element_part(e_ptr->get_local_number());

            // This should be separated out in the future, but for now provide the workspace
            e_ptr->set_workspace(&element_workspace);

            // Start calculation block
            e_ptr->start_calculations(elem_solution_vector);

            e_ptr->get_formulation().calculate(multi_field_calc);

            // End calculation block
            e_ptr->end_calculations();

            for (auto& output_algorithm : output_algorithms)
            {
                output_algorithm->add_element_contribution(*e_ptr, multi_field_calc);
            }
        }

        for (auto& output_algorithm : output_algorithms)
        {
            output_algorithm->post_loop_step();
        }

        writer_.close();
        ascii_writer_.close();
    }
}