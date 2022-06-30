#include "../UnitTestHelpers.h"
#include "bvp/bvp/FractureBVP.h"
#include "constitutive/enrichment/deformation/bonded/BondedInterfaceProperties.h"
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
      "bonded",
      constitutive::BondedInterfaceProperties::create(1e20));

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Define boundary value problem

    // A vector of { 0, 1, ..., 10 }
    std::vector<double> p(11);
    std::iota(p.begin(), p.end(), 0);
    // A vector of { 0, 1, ..., 15 }
    std::vector<double> z(16);
    std::iota(z.begin(), z.end(), 0);

    auto bvp = FractureBVP::get();
    BVPData bvp_data;
    bvp.set_input(mesh::make_linear_grid(p, p, z),
                  get_test_name(__FILE__),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);
    bvp.set_field_models(finite_element::BasalFieldSelection(),
                         "InfinitesimalStrain",
                         "isotropic",
                         bvp_data);
    add_unit_tension_z_bcs(bvp_data);

    // Insert a spherical enrichment surface defined by the radius and center below.
    const double sphere_radius1             = 4.0;
    const simple_math::Point sphere_center1 = { 5.0, 5.0, 5.0 };
    const double sphere_radius2             = 4.0;
    const simple_math::Point sphere_center2 = { 5.0, 5.0, 10.0 };
    bvp.insert_precracks_by_function_with_boolean(
      bvp_data,
      { { [sphere_radius = sphere_radius1,
           sphere_center = sphere_center1](const simple_math::Point& p) -> const double {
             const double signed_dist = simple_math::Vector3D(p[0] - sphere_center[0],
                                                              p[1] - sphere_center[1],
                                                              p[2] - sphere_center[2])
                                          .norm() -
                                        sphere_radius;
             return signed_dist;
         },
          surf_props,
          constitutive::FieldPropertiesCollectionRef() },
        { [sphere_radius = sphere_radius2,
           sphere_center = sphere_center2](const simple_math::Point& p) -> const double {
             const double signed_dist = simple_math::Vector3D(p[0] - sphere_center[0],
                                                              p[1] - sphere_center[1],
                                                              p[2] - sphere_center[2])
                                          .norm() -
                                        sphere_radius;
             return signed_dist;
         },
          surf_props,
          constitutive::FieldPropertiesCollectionRef() } });

    bvp.solve_time_steps(bvp_data);

    return 0;
}