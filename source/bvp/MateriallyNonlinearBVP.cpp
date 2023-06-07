#include "bvp/bvp/MateriallyNonlinearBVP.h"

#include "constitutive/enrichment/direction_criteria/MaxPrincipalStressDirectionCriterion.h"
#include "constitutive/enrichment/growth_criteria/AlwaysGrowGrowthCriterion.h"
#include "finite_element/FiniteElement.h"
#include "finite_element/ModelSelectionHelper.h"
#include "finite_element/hhe_primitives/PrimitiveAlgorithms.h"
#include "mesh/algorithms/Search.h"
#include "sparse_linear_algebra/LinAlgFactory.h"
#include "sparse_linear_algebra/NonLinearSolve_Newton.h"
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
    if (delta_solution != nullptr)
    {
        non_null_vector_ptrs.push_back(delta_solution.get());
    }
    if (solution != nullptr)
    {
        non_null_vector_ptrs.push_back(solution.get());
    }
    if (solution_prev != nullptr)
    {
        non_null_vector_ptrs.push_back(solution_prev.get());
    }
    if (eigenstrain_force != nullptr)
    {
        non_null_vector_ptrs.push_back(eigenstrain_force.get());
    }
    if (eigenstrain_force_prev != nullptr)
    {
        non_null_vector_ptrs.push_back(eigenstrain_force_prev.get());
    }
    if (rhs != nullptr)
    {
        non_null_vector_ptrs.push_back(rhs.get());
    }
    if (rhs_prev != nullptr)
    {
        non_null_vector_ptrs.push_back(rhs_prev.get());
    }
    if (residual != nullptr)
    {
        non_null_vector_ptrs.push_back(residual.get());
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

    if (!lin_alg.stiffness)
    {
        lin_alg.stiffness = sparse_linear_algebra::LinAlgFactory::get().create_matrix(dof_refs);
        lin_alg.stiffness->set_constraint_matrix(lin_alg.constraint_matrix.get());
        lin_alg.stiffness->allocate();
    }

    if (!lin_alg.delta_solution)
    {
        lin_alg.delta_solution =
          sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
        lin_alg.delta_solution->allocate();
    }
    if (!lin_alg.solution)
    {
        lin_alg.solution = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
        lin_alg.solution->allocate();
    }
    if (!lin_alg.solution_prev)
    {
        lin_alg.solution_prev = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
        lin_alg.solution_prev->allocate();
    }
    if (!lin_alg.eigenstrain_force)
    {
        lin_alg.eigenstrain_force =
          sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
        lin_alg.eigenstrain_force->allocate();
    }
    if (!lin_alg.eigenstrain_force_prev)
    {
        lin_alg.eigenstrain_force_prev =
          sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
        lin_alg.eigenstrain_force_prev->allocate();
    }
    if (!lin_alg.rhs)
    {
        lin_alg.rhs = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
        lin_alg.rhs->allocate();
    }
    if (!lin_alg.rhs_prev)
    {
        lin_alg.rhs_prev = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
        lin_alg.rhs_prev->allocate();
    }
    if (!lin_alg.residual)
    {
        lin_alg.residual = sparse_linear_algebra::LinAlgFactory::get().create_vector(dof_refs);
        lin_alg.residual->allocate();
    }

    if (!lin_alg.sparse_solver)
    {
        lin_alg.sparse_solver =
          sparse_linear_algebra::LinAlgFactory::get().create_solver(bvp_data.core_data.solver_opts);
    }

    _F_END();
}

void MateriallyNonlinearBVP::apply_boundary_conditions(BVPData& bvp_data) const
{
    auto& lin_alg = static_cast<MateriallyNonlinearLinearAlgebra&>(*bvp_data.gl_lin_alg);
    lin_alg.constraint_matrix->clear();

    for (auto& bc : bvp_data.bcs())
    {
        const auto scale_factor_function_index = bc->get_scale_factor_function_index();
        auto bc_scale_factor                   = 1.0;
        if (scale_factor_function_index != -1)
        {
            auto& scale_function =
              bvp_data.core_data.time_discr_data.scale_functions[scale_factor_function_index];
            const auto current_time = bvp_data.core_data.time_discr_data.current_time;
            bc_scale_factor         = scale_function.get_scale_factor(current_time) -
                              scale_function.scale_factor_for_last_converged_substep;
        }

        // IMPORTANT NOTE: this assumes the field has # of components as dimensionality
        // This works for displacement with solid elements, but reconsider what should
        // be done to handle scalar fields, shell elements, combined fields, etc.
        // TO DO figure out a robust strategy for this
        bc->initialize(bvp_data.fe_mesh(),
                       bvp_data.node_to_element_map(),
                       bvp_data.dof_manager().get_dof_maps(),
                       utility::GlobalVariables::get_dimensionality(),
                       bc_scale_factor,
                       true, /* is_incremental */
                       *lin_alg.constraint_matrix);
    }
}

void MateriallyNonlinearBVP::insert_precracks(
  BVPData& bvp_data,
  const std::vector<SurfaceInitiationAtPoint>& initations) const
{
    init_all_data(bvp_data);

    auto& enrichment_surfaces_data = bvp_data.enrichment_surfaces_data();
    auto xfem_param                = get_xfem_params(bvp_data);
    auto e_workspace               = std::make_unique<finite_element::ElementWorkspace>();

    auto& selection_helper = finite_element::ModelSelectionHelper::get();
    for (auto& initation : initations)
    {
        //////////////////////////////////////////////////////////////////////////////////////
        // Create and initiate the surface

        auto surface_model =
          selection_helper.create_enrichment_surface_model("InfinitesimalStrain",
                                                           initation.properties,
                                                           initation.growth_criterion,
                                                           initation.direction_criterion);

        auto find_result =
          mesh::find_element_with_point(bvp_data.mesh(), initation.initiation_point);
        if (!find_result.first)
        {
            throw _UTIL_MSG_EXCEPT("No element was found that contains the point " +
                                   simple_math::to_string(initation.initiation_point));
        }
        auto initation_fe = bvp_data.fe_mesh().get_finite_element(find_result.second);

        auto& direction_criteria             = initation.direction_criterion;
        simple_math::VoigtVector zero_stress = simple_math::VoigtVector::Zero(6);
        auto global_surface_id =
          static_cast<xfem::EnrichmentSurfaceID>(enrichment_surfaces_data.size() - 1);
        auto crack_normal =
          direction_criteria.get().eval(global_surface_id, initation.initiation_point, zero_stress);

        enrichment_surfaces_data.push_back(
          xfem::initiate_surface(xfem_param,
                                 surface_model.enrichment_model,
                                 surface_model.pointwise_enrichment_model_data,
                                 initation_fe,
                                 0,
                                 initation.initiation_point,
                                 crack_normal));
        xfem::update_surface_ids(xfem_param, enrichment_surfaces_data);

        xfem::write_level_sets_to_vtk(enrichment_surfaces_data,
                                      xfem_param,
                                      "./",
                                      enrichment_surfaces_data.back().global_surface_id, // step #
                                      0); // substep #
    }

    for (auto& surface_data : enrichment_surfaces_data)
    {
        //////////////////////////////////////////////////////////////////////////////////////
        // Grow the new surface until it cannot grow any further

        int growth_step = 1;
        while (!surface_data.enrichments_at_front.empty())
        {
            std::vector<xfem::EnrichmentIndexInSurface> enrichment_indices_to_grow;
            for (auto& index : surface_data.enrichments_at_front)
            {
                auto enrichment_ref = surface_data.enrichment_refs[index];
                auto front_element  = xfem::get_element(xfem_param, enrichment_ref);
                auto& enrichment    = xfem::get_enrichment(xfem_param, enrichment_ref);
                auto& criterion     = enrichment.get_growth_criterion();

                if (criterion.eval(surface_data.global_surface_id,
                                   mesh::calc_mean_coordinate(front_element),
                                   simple_math::Vector3D::Zero()))
                {
                    enrichment_indices_to_grow.push_back(index);
                }
            }

            // If no enrichments can grow due to the propagation criteria, then stop
            if (enrichment_indices_to_grow.empty())
            {
                break;
            }

            for (auto& index : enrichment_indices_to_grow)
            {
                auto search_results = xfem::search_ahead_of_front(surface_data, xfem_param, index);

                auto new_growth_occurred = xfem::grow_surface_front_into_all_leaves(surface_data,
                                                                                    xfem_param,
                                                                                    index,
                                                                                    search_results,
                                                                                    *e_workspace);
            }

            xfem::write_level_sets_to_vtk(
              enrichment_surfaces_data,
              xfem_param,
              "./",
              enrichment_surfaces_data.back().global_surface_id, // step #
              static_cast<int>(growth_step));                    // substep #

            growth_step++;
        }
    }
}

void MateriallyNonlinearBVP::renumber_dofs(BVPData& bvp_data) const
{
    if (bvp_data.dof_manager().requires_renumber(bvp_data.fe_mesh()))
    {
        auto new_to_old_dof_clone_map =
          bvp_data.dof_manager().renumber(bvp_data.fe_mesh(),
                                          bvp_data.adjacency_list(),
                                          bvp_data.element_scatter_map());

        auto& lin_alg = static_cast<MateriallyNonlinearLinearAlgebra&>(*bvp_data.gl_lin_alg);
        for (auto& v : lin_alg.vectors())
        {
            v->update_for_renumbering(new_to_old_dof_clone_map);
        }
        for (auto& m : lin_alg.matrices())
        {
            m->update_for_renumbering(new_to_old_dof_clone_map);
        }

        // Update the local copy of the solution after renumbering in case element parts are needed
        // during an assembly
        lin_alg.solution->update_local_array(bvp_data.solution());
        lin_alg.residual->update_local_array(bvp_data.core_data.residual);
    }
}

void MateriallyNonlinearBVP::insert_precracks_by_function_outside_previous(
  BVPData& bvp_data,
  const std::vector<SurfaceInsertionByLevelSet>& insertions) const
{
    init_all_data(bvp_data);

    auto& enrichment_surfaces_data = bvp_data.enrichment_surfaces_data();
    auto xfem_param                = get_xfem_params(bvp_data);

    auto always_grow_ref =
      constitutive::ensure_growth_criterion("always_grow",
                                            constitutive::AlwaysGrowGrowthCriterion());
    auto max_prin_dir_ref = constitutive::ensure_direction_criterion(
      "max_prin_stress",
      constitutive::MaxPrincipalStressDirectionCriterion());

    auto& selection_helper = finite_element::ModelSelectionHelper::get();
    for (auto i : utility::indices(insertions))
    {
        auto& insertion = insertions[i];

        auto surface_model = selection_helper.create_enrichment_surface_model("InfinitesimalStrain",
                                                                              insertion.properties,
                                                                              always_grow_ref,
                                                                              max_prin_dir_ref);

        enrichment_surfaces_data.push_back(xfem::insert_surface_based_on_level_set(
          xfem_param,
          surface_model.enrichment_model,
          surface_model.pointwise_enrichment_model_data,
          insertion.level_set_function,
          [&insertions, &i](const finite_element::FiniteElementConstAccessor& curr_fe,
                            const finite_element::LeafIndex& leaf_index,
                            const simple_math::Point& p) -> bool
          {
              // If the element is not enriched, then ensure it on the side of all of the level-sets
              if (curr_fe.layout.get_n_enrichments() == 0)
              {
                  for (auto j : utility::range(i))
                  {
                      if (insertions[j].level_set_function(p) <= 0.0)
                      {
                          return false;
                      }
                  }
                  return true;
              }
              else
              {
                  return finite_element::does_all_ancestors_match_field_sign(
                    curr_fe,
                    leaf_index,
                    finite_element::FieldSign::Pos);
              }
          },
          false,
          insertion.properties_for_fields_on_neg_side));
        xfem::update_surface_ids(xfem_param, enrichment_surfaces_data);

        xfem::write_level_sets_to_vtk(enrichment_surfaces_data,
                                      xfem_param,
                                      "./",
                                      enrichment_surfaces_data.back().global_surface_id, // step #
                                      0); // substep #
    }
}

SolveTimeStepResult MateriallyNonlinearBVP::solve_next_time_step(BVPData& bvp_data) const
{
    init_all_data(bvp_data);

    renumber_dofs(bvp_data);

    auto& time_discr_data = bvp_data.core_data.time_discr_data;

    auto current_time_step = time_discr_data.last_converged_time_step + 1;
    double target_time     = time_discr_data.target_times[current_time_step];

    std::cout << "Time step " << current_time_step << " for target time " << target_time
              << std::endl;

    // Create a directory for the current time step
    auto ts_dir = bvp_data.output_dir() / ("ts_" + std::to_string(current_time_step));
    if (!fs::exists(ts_dir))
    {
        fs::create_directory(ts_dir);
    }

    auto log_stage_name = bvp_data.bvp_name() + "_ts" + std::to_string(current_time_step);
    sparse_linear_algebra::LinAlgFactory::get().start_log_stage(log_stage_name);

    auto& lin_alg = static_cast<MateriallyNonlinearLinearAlgebra&>(*bvp_data.gl_lin_alg);

    // TO DO reevaluate how often the local array is updated and see if the tangent and residual can
    //      calculated in the same assembler execute call.
    auto& assembler = bvp_data.assembler();
    auto& dof_maps  = bvp_data.dof_manager().get_dof_maps();
    std::vector<double> local_solution(dof_maps.n_dofs_on_rank(), 0.0);
    std::vector<double> local_residual(dof_maps.n_dofs_on_rank(), 0.0);
    std::function<void(const sparse_linear_algebra::Vector&, sparse_linear_algebra::SparseMatrix&)>
      tangent_matrix_functor =
        [&local_solution, &assembler, &dof_maps](const sparse_linear_algebra::Vector& solution,
                                                 sparse_linear_algebra::SparseMatrix& A)
    {
        solution.update_local_array(local_solution);
        assembler.queue_assembly(finite_element::FEMatricesAndVectors::Flags::LHSMatrix, &A);
        assembler.execute(local_solution, dof_maps);
    };

    int substep_index = 0;
    int iteration_index = 0;
    std::function<void(const sparse_linear_algebra::Vector&, sparse_linear_algebra::Vector&)>
      residual_functor = [&local_solution,
                          &local_residual,
                          &assembler,
                          &dof_maps,
                          &bvp_data,
                          &ts_dir,
                          &substep_index,
                          &iteration_index](const sparse_linear_algebra::Vector& solution,
                                          sparse_linear_algebra::Vector& residual)
    {
        solution.update_local_array(local_solution);

        // Note: this calculates the "internal force vector"
        assembler.queue_assembly(finite_element::FEMatricesAndVectors::Flags::CompleteRHSVector,
                                 &residual);
        assembler.execute(local_solution, dof_maps);
        // TO DO implement residual = "internal force vector" (Ax) - "external force vector" (b)
        // still missing the external force part

        residual.update_local_array(local_residual);

        // TO DO user should be able to decide if iterations should be
        {
            // Save the BVP data with the current iteration solution
            bvp_data.core_data.solution.swap(local_solution);
            bvp_data.core_data.residual.swap(local_residual);

            auto substep_dir_name   = "substep_" + std::to_string(substep_index);
            auto iteration_dir_name = "it_" + std::to_string(iteration_index++);
            save_bvp_data(bvp_data,
                          (ts_dir / substep_dir_name / iteration_dir_name / "iteration_state")
                            .generic_string());

            /* Expensive debugging!
            save_bvp_data_json(bvp_data,
                               (ts_dir / substep_dir_name / iteration_dir_name / "iteration_state")
                                 .generic_string());
            */

            bvp_data.core_data.solution.swap(local_solution);
            bvp_data.core_data.residual.swap(local_residual);
        }
    };

    // Note: that these are references
    // For the first step, the value of this does not matter. For subsequent steps, it will be
    // delta_time corresponding to the last converged substep (which could be small if convergence
    // was not achieved and a cutback took place).
    double& delta_time = time_discr_data.current_delta_time;
    // The current time will be the time of the last converged step at this point. It is updated
    // inside the while loop below.
    double& current_time = time_discr_data.current_time;

    // Make sure we do not overshoot the next target time
    if ((target_time - current_time) < delta_time)
    {
        delta_time = target_time - current_time;
    }
    _UTIL_ASSERT(delta_time > 0,
                 "Somehow the logic resulted in a time increment of 0, which would result in an "
                 "infinite loop.")

    std::vector<SubstepResult> substep_results;
    while (time_discr_data.last_converged_time != target_time)
    {
        // Figure out the time that will be used for this substep
        current_time = time_discr_data.last_converged_time + delta_time;
        if (current_time > target_time)
        {
            current_time = target_time;
        }

        std::cout << "- Attempting substep from time " << time_discr_data.last_converged_time
                  << " to " << current_time << std::endl;

        // Initialize boundary conditions
        apply_boundary_conditions(bvp_data);

        // try to solve nonlinear problem
        bool substep_converged = false;
        auto status =
          sparse_linear_algebra::nonlinear_solve(bvp_data.core_data.nonlinear_solver_opts,
                                                 lin_alg.solution_prev.get(),
                                                 *lin_alg.sparse_solver,
                                                 lin_alg.stiffness.get(),
                                                 lin_alg.constraint_matrix.get(),
                                                 lin_alg.rhs.get(),
                                                 lin_alg.delta_solution.get(),
                                                 lin_alg.residual.get(),
                                                 tangent_matrix_functor,
                                                 residual_functor);
        substep_converged = status.integral_stats["converged"] == 1;
        substep_index++;
        iteration_index = 0;

        if (substep_converged)
        {
            // Update the scale factors in the scale functions so we remember the last value for the
            // calculation of the incremental scale factors.
            for (auto& scale_function : bvp_data.core_data.time_discr_data.scale_functions)
            {
                scale_function.scale_factor_for_last_converged_substep =
                  scale_function.get_scale_factor(current_time);
            }

            // Update time data
            time_discr_data.last_converged_time = current_time;

            *lin_alg.solution_prev = *lin_alg.solution;
            *lin_alg.solution      = *lin_alg.delta_solution + *lin_alg.solution_prev;
            lin_alg.solution->update_local_array(bvp_data.solution());

            /// update the internal state variables (overwritting the
            /// internal state variables for the last converged time step)
            assembler.queue_operation(
              [](finite_element::FiniteElement::ConstParams params,
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

        substep_results.push_back({ substep_converged, current_time, status });
    }

    time_discr_data.last_converged_time_step++;

    sparse_linear_algebra::LinAlgFactory::get().end_current_log_stage();

    /*
    if (BVPGlobalOptions::get().save_global_matrices)
    {
        gl_stiffness.save((ts_dir / "matrix_post_constraints.bin").generic_string());
        gl_rhs.save((ts_dir / "rhs_post_constraints.bin").generic_string());
        gl_displacement.save((ts_dir / "solution.bin").generic_string());
    }
    */

    save_bvp_data(bvp_data, (ts_dir / "converged").generic_string());

    return { true, /* converged */
             time_discr_data.current_time ==
               time_discr_data.target_times.back() /* done with time steps */,
             substep_results };
}

xfem::Params MateriallyNonlinearBVP::get_xfem_params(BVPData& bvp_data) const
{
    // Create the callback for the fea_crack module to get the element solution
    // vector when evaluating the direction criteria.
    const auto& solution = bvp_data.solution();
    const auto& dof_map  = bvp_data.dof_manager().get_dof_maps();
    xfem::GetElementSolutionCallback get_element_solution_cb =
      [&solution, &dof_map](
        const finite_element::FiniteElementConstAccessor& fe_accessor) -> simple_math::GeneralVector
    {
        return assembler::extract_element_solution_vector(solution,
                                                          dof_map,
                                                          fe_accessor.element.local_index());
    };

    return xfem::make_params(bvp_data.fe_mesh(),
                             bvp_data.adjacency_list(),
                             get_element_solution_cb,
                             bvp_data.dof_manager());
}

} // end namespace bvp
