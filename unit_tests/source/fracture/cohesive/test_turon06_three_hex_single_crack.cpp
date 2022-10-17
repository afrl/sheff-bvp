#include "../../UnitTestHelpers.h"
#include "bvp/bvp/MateriallyNonlinearBVP.h"
#include "constitutive/enrichment/deformation/turon_06_cz/InfStrainTuron06EnrichmentModel.h"
#include "constitutive/enrichment/deformation/turon_06_cz/Turon06EnrichmentProperties.h"
#include "constitutive/enrichment/direction_criteria/PositionFunctionDirectionCriterion.h"
#include "constitutive/enrichment/growth_criteria/AlwaysGrowGrowthCriterion.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "mesh/algorithms/Creation.h"

using namespace bvp;

int main(int argc, char** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::get().init(
      sparse_linear_algebra::LinAlgFactory::LinAlgBackend::AUTO,
      {});

    bvp::BVPGlobalOptions::get().save_global_matrices = true;

    // 1 element problem for which we know the closed solution
    {
        //////////////////////////////////////////////////////////////////////////////////////////////
        // Properties

        constitutive::add_field_properties(
          "isotropic",
          constitutive::LinearElasticProperties::create_isotropic_E_nu(3.0,
                                                                       0.0 /*force 1d solution*/));

        auto surf_props = constitutive::add_enrichment_properties(
          "cohesive",
          constitutive::Turon06EnrichmentProperties::create(9.0 /* J/m^2 */,
                                                            9.0 /* J/m^2 */,
                                                            1.0 /* eta */,
                                                            3.0 /* Pa */,
                                                            3.0 /* Pa */,
                                                            3.0 /* Pa/m */));


        auto direction_crit = constitutive::add_direction_criterion(
          "direction_crit",
          constitutive::PositionFunctionDirectionCriterion(
            [](const int& crack_index, const simple_math::Point& p) -> simple_math::Vector3D
            { return simple_math::Vector3D(1, 0, 0); }));

        auto growth_crit =
          constitutive::add_growth_criterion("growth_crit",
                                             constitutive::AlwaysGrowGrowthCriterion());

        //////////////////////////////////////////////////////////////////////////////////////////////
        // Define boundary value problem

        auto bvp = MateriallyNonlinearBVP::get();
        BVPData bvp_data;
        bvp.set_input(mesh::make_linear_grid({ 0, 2 }, { 0, 1 }, { 0, 1 }),
                      get_test_name(__FILE__),
                      (get_output_directory(__FILE__) / "output").string(),
                      bvp_data);
        bvp.set_field_models(finite_element::BasalFieldSelection(),
                             "InfinitesimalStrain",
                             "isotropic",
                             bvp_data);
        add_symmetry_about_0_planes_bcs(bvp_data);
        bvp_data.add_bc(DirichletBC("pos_x", { 0, 4.0 /* U */ }));

        bvp.insert_precracks(bvp_data,
                             { { { 1.0, 0.5, 0.5 }, // point
                                 surf_props,
                                 growth_crit,
                                 direction_crit } });

        bvp_data.core_data.nonlinear_solver_opts.rel_residual_norm_tol = 1e-6;
        auto solve_result = bvp.solve_next_time_step(bvp_data);
    }

    return 0;

    // 3 element problem
    {
        //////////////////////////////////////////////////////////////////////////////////////////////
        // Properties

        constitutive::add_field_properties(
          "isotropic",
          constitutive::LinearElasticProperties::create_isotropic_E_nu(1e9, 0.3));

        auto surf_props = constitutive::add_enrichment_properties(
          "cohesive",
          constitutive::Turon06EnrichmentProperties::create(177 /* J/m^2 */,
                                                            496 /* J/m^2 */,
                                                            2.25,
                                                            77e6 /* MPa */,
                                                            100e6 /* MPa */,
                                                            1e20 /* MPa */));

        auto direction_crit = constitutive::add_direction_criterion(
          "direction_crit",
          constitutive::PositionFunctionDirectionCriterion(
            [](const int& crack_index, const simple_math::Point& p) -> simple_math::Vector3D
            { return simple_math::Vector3D(1, 0, 0); }));

        auto growth_crit =
          constitutive::add_growth_criterion("growth_crit",
                                             constitutive::AlwaysGrowGrowthCriterion());

        //////////////////////////////////////////////////////////////////////////////////////////////
        // Define boundary value problem

        auto bvp = MateriallyNonlinearBVP::get();
        BVPData bvp_data;
        bvp.set_input(mesh::make_linear_grid({ 0, 1, 2, 3 }, { 0, 1 }, { 0, 1 }),
                      get_test_name(__FILE__),
                      (get_output_directory(__FILE__) / "output").string(),
                      bvp_data);
        bvp.set_field_models(finite_element::BasalFieldSelection(),
                             "InfinitesimalStrain",
                             "isotropic",
                             bvp_data);
        add_unit_tension_x_bcs(bvp_data);

        bvp.insert_precracks(bvp_data,
                             { { { 1.525, 0.2, 0.2 }, // point
                                 surf_props,
                                 growth_crit,
                                 direction_crit } });

        bvp_data.core_data.nonlinear_solver_opts.rel_residual_norm_tol = 1e-6;
        bvp.solve_time_steps(bvp_data);
    }

    return 0;
}