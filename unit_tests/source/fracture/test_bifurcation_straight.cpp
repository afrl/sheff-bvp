#include "../UnitTestHelpers.h"
#include "bvp/bvp/FractureBVP.h"
#include "constitutive/enrichment/deformation/traction_free/TractionFreeInterfaceProperties.h"
#include "constitutive/enrichment/direction_criteria/PositionFunctionDirectionCriterion.h"
#include "constitutive/enrichment/growth_criteria/PositionFunctionGrowthCriterion.h"
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

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    const double E  = 1e9;
    const double nu = 0.3;
    constitutive::add_field_properties(
      "isotropic",
      constitutive::LinearElasticProperties::create_isotropic_E_nu(E, nu));

    auto surf_props = constitutive::add_enrichment_properties(
      "traction_free",
      constitutive::TractionFreeInterfaceProperties::create());

    auto crack_0_dir_crit = constitutive::add_direction_criterion(
      "crack_0_dir_crit",
      constitutive::PositionFunctionDirectionCriterion(
        [](const int& crack_index, const simple_math::Point& p) -> simple_math::Vector3D {
            return simple_math::Vector3D(0, 0, 1);
        }));

    auto crack_1_dir_crit = constitutive::add_direction_criterion(
      "crack_1_dir_crit",
      constitutive::PositionFunctionDirectionCriterion(
        [](const int& crack_index, const simple_math::Point& p) -> simple_math::Vector3D {
            return simple_math::Vector3D(1, 0, 0);
        }));

    int phase        = 0;
    auto growth_crit = constitutive::add_growth_criterion(
      "growth_crit",
      constitutive::PositionFunctionGrowthCriterion(
        [&phase](const int& crack_index, const simple_math::Point& p) -> bool {
            switch (phase)
            {
                case 0: return (p[1] < 4 && crack_index == 0) || (p[1] > 6 && crack_index == 1);
                case 1: return crack_index == 0;
                case 2: return crack_index == 1;
                default: return false;
            }
        }));

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Define boundary value problem

    // A vector of { 0, 1, ..., 10 }
    std::vector<double> p(11);
    std::iota(p.begin(), p.end(), 0);

    auto bvp = FractureBVP::get();
    BVPData bvp_data;
    bvp.set_input(mesh::make_linear_grid(p, p, p),
                  get_test_name(__FILE__),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);
    bvp.set_field_models(finite_element::BasalFieldSelection(),
                         "InfinitesimalStrain",
                         "isotropic",
                         bvp_data);
    add_unit_tension_x_z_bcs(bvp_data);

    bvp.insert_precracks(bvp_data,
                         { { { 4.5, 0.2, 4.5 }, // point
                             surf_props,
                             growth_crit,
                             crack_0_dir_crit },
                           { { 4.5, 9.8, 4.5 }, // point
                             surf_props,
                             growth_crit,
                             crack_1_dir_crit } });

    // Let the cracks grow past each other
    phase = 1;
    bvp.try_to_grow_precracks(bvp_data);
    phase = 2;
    bvp.try_to_grow_precracks(bvp_data);

    bvp.solve_time_steps(bvp_data);

	return 0;
}