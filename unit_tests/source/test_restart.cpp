#include "UnitTestHelpers.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "bvp/bvp/LinearElasticBVP.h"
#include "fea_io/common/MPIUtils.h"
#include "mesh/algorithms/UnitElements.h"
#include "sparse_linear_algebra/LinAlgFactory.h"
#include "utility/MPIUtils.h"

using namespace bvp;

int main(int argc, char** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::init(
      sparse_linear_algebra::LinAlgFactory::LinAlgBackend::AUTO,
      &argc,
      &argv);

    // Create a constitutive relation
    const double E  = 1e9;
    const double nu = 0.3;
    auto props = constitutive::add_field_properties(
      "isotropic",
      constitutive::LinearElasticProperties::create_isotropic_E_nu(E, nu));

    auto bvp = LinearElasticBVP::get();
    BVPData bvp_data;

    bvp.set_input(mesh::get_unit_linear_hex(),
                  get_test_name(__FILE__),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);
    bvp.set_field_models(finite_element::BasalFieldSelection(),
                         "InfinitesimalStrain",
                         "isotropic",
                         bvp_data);
    add_unit_tension_x_bcs(bvp_data);
    bvp.solve_time_steps(bvp_data);

    save_bvp_data_json(bvp_data, (get_output_directory(__FILE__) / "converged").generic_string());

    auto recovered_bvp_data =
     load_bvp_data((get_output_directory(__FILE__) / "output/ts_0/converged").generic_string());

    save_bvp_data_json(recovered_bvp_data,
                       (get_output_directory(__FILE__) / "recovered").generic_string());

    return 0;
}
