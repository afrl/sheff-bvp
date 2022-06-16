#include "../UnitTestHelpers.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "bvp/bvp/LinearElasticBVP.h"
#include "mesh/algorithms/UnitElements.h"

using namespace bvp;

int main(int argc, char** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::get().init(
      sparse_linear_algebra::LinAlgFactory::LinAlgBackend::AUTO,
      {});

    const double E  = 1e9;
    const double nu = 0.3;
    constitutive::add_field_properties(
      "isotropic",
      constitutive::LinearElasticProperties::create_isotropic_E_nu(E, nu));

    auto bvp = LinearElasticBVP::get();
    BVPData bvp_data;
    bvp.set_input(mesh::get_unit_element(mesh::ElementType::H3D8N),
                  get_test_name(__FILE__),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);
    bvp.set_field_models(finite_element::BasalFieldSelection(),
                         "InfinitesimalStrain",
                         "isotropic",
                         bvp_data);
    add_unit_tension_x_bcs(bvp_data);
    bvp.solve_time_steps(bvp_data);

    return 0;
}
