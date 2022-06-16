#include "../UnitTestHelpers.h"
#include "bvp/bvp/FractureBVP.h"
#include "constitutive/enrichment/deformation/traction_free/TractionFreeInterfaceProperties.h"
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

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    constitutive::add_field_properties(
      "isotropic",
      constitutive::LinearElasticProperties::create_isotropic_E_nu(1e9, 0.3));

    auto surf_props = constitutive::add_enrichment_properties(
      "traction_free",
      constitutive::TractionFreeInterfaceProperties::create());

    auto crack_0_dir_crit = constitutive::add_direction_criterion(
      "crack_0_dir_crit",
      constitutive::PositionFunctionDirectionCriterion(
        [](const int& crack_index, const simple_math::Point& p) -> simple_math::Vector3D {
            return simple_math::Vector3D(1, 0, 0);
        }));

    auto growth_crit =
      constitutive::add_growth_criterion("growth_crit", constitutive::AlwaysGrowGrowthCriterion());

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Define boundary value problem

    auto N = 9;
    std::vector<double> grid_offsets(N);
    for (auto i : utility::indices(grid_offsets))
    {
        grid_offsets[i] = 2.0 * double(i) / double(grid_offsets.size() - 1) - 1.0;
    }

    auto bvp = FractureBVP::get();
    BVPData bvp_data;
    bvp.set_input(mesh::make_linear_grid(grid_offsets, grid_offsets, { 0, 0.2, 0.4 }),
                  get_test_name(__FILE__),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);
    bvp.set_field_models(finite_element::BasalFieldSelection(),
                         "InfinitesimalStrain",
                         "isotropic",
                         bvp_data);
    add_unit_tension_x_z_bcs(bvp_data);

    bvp.insert_precracks(bvp_data,
                         { { { 0.525, 0.2, 0.2 }, // point
                             surf_props,
                             growth_crit,
                             crack_0_dir_crit } });

    bvp.solve_time_steps(bvp_data);

    return 0;
}