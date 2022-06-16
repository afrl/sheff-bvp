#include "../UnitTestHelpers.h"
#include "bvp/bvp/FractureBVP.h"
#include "constitutive/enrichment/deformation/traction_free/TractionFreeInterfaceProperties.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"

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

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Create extra node sets

    // A vector of { 0, 1, ..., 10 }
    std::vector<double> p(11);
    std::iota(p.begin(), p.end(), 0);

    auto grid = mesh::make_linear_grid(p, p, p);
    mesh::write_vtk(grid, "grid.vtk", false);

    mesh::add_to_group(grid,
                       "fixed_sphere_nodes",
                       mesh::NodeSelection([](const mesh::NodeConstAccessor& n) -> bool {
                           return (n.x() == 5.0 && n.y() == 5.0 && n.z() == 5.0) ||
                                  (n.x() == 6.0 && n.y() == 5.0 && n.z() == 5.0) ||
                                  (n.x() == 5.0 && n.y() == 6.0 && n.z() == 5.0);
                       }));

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Define boundary value problem


    auto bvp = FractureBVP::get();
    BVPData bvp_data;
    bvp.set_input(std::move(grid),
                  get_test_name(__FILE__),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);
    bvp.set_field_models(finite_element::BasalFieldSelection(),
                         "InfinitesimalStrain",
                         "isotropic",
                         bvp_data);

    add_unit_tension_z_bcs(bvp_data);
    // Add an extra boundary condition to move the sphere in the center upwards and prevent rigid
    // body rotation.
    bvp_data.add_bc(DirichletBC("fixed_sphere_nodes", { 0, 0.0, BCApplicationType::NODE }));
    bvp_data.add_bc(DirichletBC("fixed_sphere_nodes", { 1, 0.0, BCApplicationType::NODE }));
    bvp_data.add_bc(DirichletBC("fixed_sphere_nodes", { 2, 0.5, BCApplicationType::NODE }));

    // 1) Insert a spherical enrichment surface defined by the radius and center below.
    // 2) Insert a curved enrichment surface defined by the radius and center below that
    //    is initiated outside the sphere.
    const double sphere_1_radius             = 3.7;
    const simple_math::Point sphere_1_center = { 5.0, 5.0, 5.0 };
    const double sphere_2_radius             = 8.0;
    simple_math::Point sphere_2_center       = { 5, 5, 12 };
    constitutive::TractionFreeInterfaceProperties enrichment_properties;
    bvp.insert_precracks_by_function_outside_previous(
      bvp_data,
      { { // surface 0
          [sphere_1_radius, sphere_1_center](const simple_math::Point& p) -> const double {
              const double signed_dist = simple_math::Vector3D(p[0] - sphere_1_center[0],
                                                               p[1] - sphere_1_center[1],
                                                               p[2] - sphere_1_center[2])
                                           .norm() -
                                         sphere_1_radius;
              return signed_dist;
          },
          surf_props,
          constitutive::FieldPropertiesCollectionRef() },
        { [sphere_2_radius, sphere_2_center](const simple_math::Point& p) -> const double {
             const double signed_dist = simple_math::Vector3D(p[0] - sphere_2_center[0],
                                                              p[1] - sphere_2_center[1],
                                                              p[2] - sphere_2_center[2])
                                          .norm() -
                                        sphere_2_radius;
             return signed_dist;
         },
          surf_props,
          constitutive::FieldPropertiesCollectionRef() } });

    bvp.solve_time_steps(bvp_data);

    return 0;
}