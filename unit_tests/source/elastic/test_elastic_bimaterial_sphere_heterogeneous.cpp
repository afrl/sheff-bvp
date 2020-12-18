#include "../UnitTestHelpers.h"
#include "bvp/bvp/LinearElasticBVP.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "mesh/algorithms/IO.h"

using namespace bvp;

int main(int argc, char** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::init(
      sparse_linear_algebra::LinAlgFactory::LinAlgBackend::AUTO,
      &argc,
      &argv);

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    constitutive::add_field_properties(
      "matrix",
      constitutive::LinearElasticProperties::create_isotropic_E_nu(1e9, 0.3));

    constitutive::add_field_properties(
      "inclusion",
      constitutive::LinearElasticProperties::create_isotropic_E_nu(1e9, 0.3));

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Define boundary value problem

    auto bvp = LinearElasticBVP::get();
    BVPData bvp_data;
    bvp.set_input(mesh::read_abaqus(get_test_data_filepath("spherical_inclusion_refined.inp")),
                  get_test_name(__FILE__),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);

    // Set properties and model for matrix elements
    bvp.set_field_models(
      finite_element::BasalFieldSelection([](const mesh::ElementConstAccessor& element,
                                             const finite_element::LeafIndex& basal_field_index) {
          return element.component_id() == 1;
      }),
      "InfinitesimalStrain",
      "matrix",
      bvp_data);

    // Set properties and model for inclusion elements
    bvp.set_field_models(
      finite_element::BasalFieldSelection([](const mesh::ElementConstAccessor& element,
                                             const finite_element::LeafIndex& basal_field_index) {
          return element.component_id() == 0;
      }),
      "InfinitesimalStrain",
      "inclusion",
      bvp_data);

    add_unit_tension_z_bcs(bvp_data);

    bvp.solve_time_steps(bvp_data);

    return 0;
}