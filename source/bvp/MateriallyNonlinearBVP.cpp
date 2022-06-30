#include "bvp/bvp/MateriallyNonlinearBVP.h"

#include "finite_element/FiniteElement.h"
#include "sparse_linear_algebra/LinAlgFactory.h"
#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/Vector.h"
#include "utility/Profile.h"
#include "utility/serialization/Serialization.h"

#include <chrono>

namespace bvp
{

LoadStepManagerOpts loadParameters;
/// parameters have to come form some input
/// The maximum number of iterations allowed for any iteration loop.
loadParameters.max_iterations=100;
  /// step reduction if convergence not achieved
loadParameters.reduction_coeff=2;
  /// step increase when convergence is achieved in a number of steps
loadParamters.increase_coeff=1.2;
  /// number of converged steps to start increasing the step
loadParameters.numberConvergedStepsForIncrease;
  /// load sequence from the input file????
loadParameters.load_sequence_init;
  ///
loadParameters.loadSequenceIndex;

NonlinearStaticLoadStepManager  loadSequence(const LoadStepManagerOpts loadParameters);
///
/// this load sequence is going to change during solution
///
loadSequence.load_sequence_init_current = loadParameters.load_sequence_init;

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

    lin_alg.solution = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.solution->allocate();
    lin_alg.solution_prev = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.solution_prev->allocate();
    lin_alg.eigenstrain_force = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.eigenstrain_force->allocate();
    lin_alg.eigenstrain_force_prev = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.eigenstrain_force_prev->allocate();
    lin_alg.rhs = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.rhs->allocate();
    lin_alg.rhs_prev = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
    lin_alg.rhs_prev->allocate();

    lin_alg.sparse_solver =
      sparse_linear_algebra::LinAlgFactory::get().create_solver(bvp_data.core_data.solver_opts);

    _F_END();
}

void MateriallyNonlinearBVP::apply_boundary_conditions(BVPData& bvp_data, float loadCoefficient) const
{
    auto& lin_alg = static_cast<MateriallyNonlinearLinearAlgebra&>(*bvp_data.gl_lin_alg);

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
                       loadCoefficient,
                       *lin_alg.constraint_matrix);
    }
}

MateriallyNonlinearBVP::solve_time_steps(BVPData& bvp_data) const
//
//
{   
    while (not_last_step_done()) {
        _F_BEGIN2(bvp_data.bvp_name());

        loadSequence.increase_loadSequenceIndex();
      
        init_all_data(bvp_data);
    
        // Write restart data
        // Create a directory for the current time step
        if (loadSequence.get_previous_step_state()) {
            bvp_data.core_data.time_discr_data.last_converged_time_step++;
            auto curr_time_step = bvp_data.core_data.time_discr_data.last_converged_time_step;
            auto ts_dir         = bvp_data.output_dir() / ("ts_" + std::to_string(curr_time_step));
            if (!fs::exists(ts_dir))
            {
                fs::create_directory(ts_dir);
            }
        
            auto log_stage_name = bvp_data.bvp_name() + "_ts" + std::to_string(curr_time_step);
        // After residual is converged, update the internal state variables (overwritting the
        // internal state variables for the last converged time step)
            sparse_linear_algebra::LinAlgFactory::get().start_log_stage(log_stage_name);
        }
        _F_BEGIN2("apply BC's");
    


        // TODO indicate if Secant or Tangent with BCUpdateType::Secant
        apply_boundary_conditions(bvp_data, loadSequence.loadCoefficient);
    
        _F_END();
    
        //////////////
        // Assembly //
        //////////////

        if (current_load_step()==1) {
            _F_BEGIN2("assembly");
            auto& assembler            = bvp_data.assembler();
            auto& lin_alg              = static_cast<MateriallyNonlinearLinearAlgebra&>(*bvp_data.gl_lin_alg);
            auto& gl_stiffness         = *lin_alg.stiffness;
            auto& gl_disp_increment    = *lin_alg.solution;
            auto& gl_displacement      = *lin_alg.solution_prev;
            auto& gl_eigenstrain_force = *lin_alg.eigenstrain_force;
            auto& gl_eigen_force_prev  = *lin_alg.eigenstrain_force_prev;
            auto& gl_rhs               = *lin_alg.rhs;
            auto& gl_rhs_prev          = *lin_alg.rhs_prev;

            gl_rhs_prev  = 0.0;
            gl_eigen_force_prev = 0.0;
            gl_displacement = 0.0;
    
        // Queue assembly of the global stiffness matrix
            assembler.queue_assembly(finite_element::FEMatricesAndVectors::Flags::LHSMatrix, &gl_stiffness);
    
        // Queue assembly of the global vector containing RHS terms pertaining
        // to eigenstrains (such as thermal strain, plastic strain, etc.)
            assembler.queue_assembly(finite_element::FEMatricesAndVectors::Flags::RHSVector,
                                 &gl_eigenstrain_force);
    
        // Execute both assemblies as one element loop
            assembler.execute(bvp_data.solution(), bvp_data.dof_manager().get_dof_maps());
    
            _F_END();
    
        // Total RHS = constraint RHS + applied loads + inelastic RHS (temperature/moisture terms)
            gl_rhs = gl_rhs + gl_eigenstrain_force;
    
        if (BVPGlobalOptions::get().save_global_matrices)
            {
                gl_stiffness.save((ts_dir / "matrix_pre_constraints.bin").generic_string());
                gl_rhs.save((ts_dir / "rhs_pre_constraints.bin").generic_string());
            }
        ///
        ///   this is a call to the NonlinearSlover that Doug develops; not working yet;
        ///   
        ///   Note: the function returns displacements, not incremental displacements!!!!
        ///
            SolverStatus sparse_linear_algebra::nonlinear_solve(
                const NonLinearSolverParameters &parameters,
                &gl_displacement,
                SparseSolver& solver, 
                &gl_stiffness, 
                &gl_rhs, 
                &gl_disp_increment,
                std::function<void(const Vector&, SparseMatrix&)> update_tangent_matrix,
                std::function<void(const Vector&, Vector&)> update_residual,
                std::function<bool(const Vector&, double&, unsigned int&)> convergence_criterion);
        
            if (state.integral_stats["converged"] ==1) {
                /// update solution vectors/// step increase when convergence is achieved in a number of steps (if previously it was decreased)
                gl_rhs_prev  = gl_rhs;
                gl_eigen_force_prev = gl_eigenstrain_force;
                gl_displacement = gl_disp_increment;
                
                sparse_linear_algebra::LinAlgFactory::get().end_current_log_stage();
        
                /// update the internal state variables (overwritting the
                /// internal state variables for the last converged time step)
                assembler.queue_operation(
                        [](finite_element::FiniteElement::ConstParams& params,
                        finite_element::FiniteElementAccessor& fe)
                        { finite_element::FiniteElement::update_converged_isvs(params, fe.layout); });
            
                // Execute the update for ISVs
                assembler.execute(/* updated solution goes here */ bvp_data.solution(),
                                    bvp_data.dof_manager().get_dof_maps());
            
                bool convergence_achieved = true;
                if (loadSequence.numberConvergedStepsForIncrease>0){
                    loadSequence.numberConvergedStepsForIncrease += 1;
                }
            }
            else {
                /// not converged; make the step smaller
                /// 
                loadSequence.add_additional_steps();
            }
        }
        else {
             _F_BEGIN2("assembly");
            auto& assembler            = bvp_data.assembler();
            auto& lin_alg              = static_cast<MateriallyNonlinearLinearAlgebra&>(*bvp_data.gl_lin_alg);
            auto& gl_stiffness         = *lin_alg.stiffness;
            auto& gl_disp_increment    = *lin_alg.solution;
            auto& gl_displacement      = *lin_alg.solution_prev;
            auto& gl_eigenstrain_force = *lin_alg.eigenstrain_force;
            auto& gl_eigen_force_prev  = *lin_alg.eigenstrain_force_prev;
            auto& gl_rhs               = *lin_alg.rhs;
            auto& gl_rhs_prev          = *lin_alg.rhs_prev;
    
        // Queue assembly of the global stiffness matrix
            assembler.queue_assembly(finite_element::FEMatricesAndVectors::Flags::LHSMatrix, &gl_stiffness);
    
        // Queue assembly of the global vector containing RHS terms pertaining
        // to eigenstrains (such as thermal strain, plastic strain, etc.)
            assembler.queue_assembly(finite_element::FEMatricesAndVectors::Flags::RHSVector,
                                 &gl_eigenstrain_force);
    
        // Execute both assemblies as one element loop
            assembler.execute(bvp_data.solution(), bvp_data.dof_manager().get_dof_maps());
    
            _F_END();
    
        // Total RHS = constraint RHS + applied loads + inelastic RHS (temperature/moisture terms)
            gl_rhs = gl_rhs + gl_eigenstrain_force - gl_eigen_force_prev;
    
        if (BVPGlobalOptions::get().save_global_matrices)
            {
                gl_stiffness.save((ts_dir / "matrix_pre_constraints.bin").generic_string());
                gl_rhs.save((ts_dir / "rhs_pre_constraints.bin").generic_string());
            }
        ///
        ///   this is a call to the NonlinearSlover that Doug develops; not working yet;
        ///   
        ///   Note: the function returns displacements, not incremental displacements!!!!
        ///
            SolverStatus sparse_linear_algebra::nonlinear_solve(
                const NonLinearSolverParameters &parameters,
                gl_displacement,
                SparseSolver& solver, 
                &gl_stiffness, 
                &gl_rhs, 
                &gl_disp_increment,
                std::function<void(const Vector&, SparseMatrix&)> update_tangent_matrix,
                std::function<void(const Vector&, Vector&)> update_residual,
                std::function<bool(const Vector&, double&, unsigned int&)> convergence_criterion);
        
            if (state.integral_stats["converged"] ==1) {
                /// update solution vectors/// step increase when convergence is achieved in a number of steps (if previously it was decreased)
                gl_rhs_prev  = gl_rhs;
                gl_eigen_force_prev = gl_eigenstrain_force;
                gl_displacement = gl_displacement + gl_disp_increment;
                
                sparse_linear_algebra::LinAlgFactory::get().end_current_log_stage();
        
                /// update the internal state variables (overwritting the
                /// internal state variables for the last converged time step)
                assembler.queue_operation(
                        [](finite_element::FiniteElement::ConstParams& params,
                        finite_element::FiniteElementAccessor& fe)
                        { finite_element::FiniteElement::update_converged_isvs(params, fe.layout); });
            
                // Execute the update for ISVs
                assembler.execute(/* updated solution goes here */ bvp_data.solution(),
                                    bvp_data.dof_manager().get_dof_maps());
            
                bool convergence_achieved = true;
                loadSequence.loadSequenceIndex += 1;
                loadSequence.number_converged_steps += 1; 
                if (loadSequence.numberConvergedStepsForIncrease>0){
                    loadSequence.numberConvergedStepsForIncrease += 1;
                }
            }
            else {
                /// not converged; make the step smaller
                /// 
                loadSequence.add_additional_steps();
            }
        }
    }
}        

    NonlinearStaticLoadStepManager::increase_loadSequenceIndex() const{
        loadSequenceIndex += 1;
        loadCoefficient += load_sequence_init_current(loadSequenceIndex);
    }
    /// 
    NonlinearStaticLoadStepManager::decrease_loadSequenceIndex() const{
        loadSequenceIndex -= 1;
        loadCoefficient -= load_sequence_init_current(loadSequenceIndex);
    }
    /// returns true if it is not over the last step
    NonlinearStaticLoadStepManage::not_last_step_done() {
        if (loadSequenceIndex <= sizeof(load_sequence_init_current) ){
            return true;
        }
        else {
            return false;
        }
    }
    /// returns max_iterations  
    NonlinearStaticLoadStepManage::get_max_iterations(){
        return LSopts_.max_iterations;
    }
    /// step reduction if convergence not achieved
    NonlinearStaticLoadStepManage::get_reduction_coeff() {
        return LSopts_.reduction_coeff;
    }
    /// returns the step increase coefficient 
    NonlinearStaticLoadStepManage::get_increase_coeff(){
        return LSopts_.increase_coeff;
    }
    /// returns the number of converged steps to start increasing the step
    NonlinearStaticLoadStepManage::get_numberConvergedStepsForIncrease(){
        return numberConsecutiveConvergedSteps;
    }
    /// returns state of the previous step - converged or not and changes the state
    NonlinearStaticLoadStepManage::get_previous_step_state(){
        return previous_step_converged;
    }
    NonlinearStaticLoadStepManage::add_additional_steps(){
        if (loadSequenceIndex > 0 ) {
            float stepSize = load_sequence_init_current(loadSequenceIndex) - load_sequence_init_current(loadSequenceIndex-1);
        }
        else {
           float stepSize = load_sequence_init_current(0); 
        }
        stepSize = stepSize/LSopts_.reduction_coeff;
        for (int i=0, i < LSopts_.reduction_coeff, i++) {
                int insert_index = begin()+loadSequenceIndex + i - 1;
                bvp_data.core_data.load_step.insert(insert_index,stepSize;
        }
        numberConvergedStepsForIncrease= 0;
    }

} // end namespace bvp
