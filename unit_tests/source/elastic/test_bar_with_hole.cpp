#include "../UnitTestHelpers.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "bvp/bvp/LinearElasticBVP.h"

using namespace bvp;

int main(int argc, char ** argv)
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

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Define boundary value problem

    auto bvp = LinearElasticBVP::get();
    BVPData bvp_data;
    bvp.set_input("",
                  get_test_data_filepath("bar_with_hole"),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);
    bvp.set_field_models(finite_element::BasalFieldSelection(),
                         "InfinitesimalStrain",
                         "isotropic",
                         bvp_data);
    
    add_extent_node_groups(bvp_data.mesh());
    bvp_data.add_bc(DirichletBC("neg_x", { 0, 0 }));
    bvp_data.add_bc(DirichletBC("neg_x", { 1, 0 }));
    bvp_data.add_bc(DirichletBC("neg_x", { 2, 0 }));
    bvp_data.add_bc(DirichletBC("pos_x", { 0, 10 }));

    bvp.solve_time_steps(bvp_data);

    return 0;
}