#include "../UnitTestHelpers.h"
#include "bvp/bvp/FractureBVP.h"
#include "constitutive/enrichment/deformation/bonded/BondedInterfaceProperties.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "mesh/algorithms/Creation.h"
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

    // Write matrices (and vectors) to binary files to examine.
    BVPGlobalOptions::get().save_global_matrices = true;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    auto matrix = constitutive::add_field_properties(
      "matrix",
      constitutive::LinearElasticProperties::create_isotropic_E_nu(1e9, 0.3));

    auto surf_props = constitutive::add_enrichment_properties(
      "bonded",
      constitutive::BondedInterfaceProperties::create(1e20));

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Define boundary value problems

    auto bvp = FractureBVP::get();
    for (auto& circle_radius : { 2.0, 2.1, 2.2, 2.3, 2.4, 2.5 })
    {
        std::string output_dir =
          (get_output_directory(__FILE__) / ("output_r_" + std::to_string(circle_radius)))
            .generic_string();

        BVPData bvp_data;
        bvp.set_input(mesh::make_linear_grid({ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 },
                                             { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 },
                                             { 0, 1 }),
                      get_test_name(__FILE__),
                      output_dir,
                      bvp_data);
        bvp.set_field_models(finite_element::BasalFieldSelection(),
                             "InfinitesimalStrain",
                             "matrix",
                             bvp_data);

        // Add simple boundary conditions to remove rigid body rotation/translation
        mesh::add_to_group(bvp_data.mesh(),
                           "fixed_nodes",
                           mesh::NodeSelection([](const mesh::NodeConstAccessor& n) -> bool {
                               return n.point() == simple_math::Point({ 0.0, 0.0, 0.0 }) ||
                                      n.point() == simple_math::Point({ 9.0, 0.0, 0.0 }) ||
                                      n.point() == simple_math::Point({ 9.0, 9.0, 0.0 });
                           }));
        assert(bvp_data.mesh().nodes().group("fixed_nodes").node_indices.size() == 3);
        bvp_data.add_bc(DirichletBC("fixed_nodes", { 0, 0, BCApplicationType::NODE }));
        bvp_data.add_bc(DirichletBC("fixed_nodes", { 1, 0, BCApplicationType::NODE }));
        bvp_data.add_bc(DirichletBC("fixed_nodes", { 2, 0, BCApplicationType::NODE }));

        const simple_math::Point circle_center = { 4.5, 4.5, 0 };
        bvp.insert_precracks_by_function_outside_previous(
          bvp_data,
          { { // surface 0
              [circle_radius, circle_center](const simple_math::Point& p) -> const double {
                  const double signed_dist =
                    simple_math::Vector3D(p[0] - circle_center[0], p[1] - circle_center[1], 0.0)
                      .norm() -
                    circle_radius;
                  return signed_dist;
              },
              surf_props,
              constitutive::FieldPropertiesCollectionRef() } });

        bvp.solve_time_steps(bvp_data);

        xfem::write_dofs_to_vtk(bvp_data.enrichment_surfaces_data(),
                                xfem::make_params(bvp_data.fe_mesh(),
                                                  bvp_data.adjacency_list(),
                                                  xfem::GetElementSolutionCallback(),
                                                  bvp_data.dof_manager()),
                                output_dir,
                                0);
        xfem::write_level_sets_to_vtk(bvp_data.enrichment_surfaces_data(),
                                      xfem::make_params(bvp_data.fe_mesh(),
                                                        bvp_data.adjacency_list(),
                                                        xfem::GetElementSolutionCallback(),
                                                        bvp_data.dof_manager()),
                                      output_dir,
                                      0);
        mesh::write_vtk(bvp_data.mesh(), output_dir + "/mesh.vtk");
    }

    {
        std::string output_dir =
          (get_output_directory(__FILE__) / ("output_no_surface")).generic_string();

        BVPData bvp_data;
        bvp.set_input(mesh::make_linear_grid({ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 },
                                             { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 },
                                             { 0, 1 }),
                      get_test_name(__FILE__),
                      output_dir,
                      bvp_data);
        bvp.set_field_models(finite_element::BasalFieldSelection(),
                             "InfinitesimalStrain",
                             "matrix",
                             bvp_data);

        // Add simple boundary conditions to remove rigid body rotation/translation
        mesh::add_to_group(bvp_data.mesh(),
                           "fixed_nodes",
                           mesh::NodeSelection([](const mesh::NodeConstAccessor& n) -> bool {
                               return n.point() == simple_math::Point({ 0.0, 0.0, 0.0 }) ||
                                      n.point() == simple_math::Point({ 9.0, 0.0, 0.0 }) ||
                                      n.point() == simple_math::Point({ 9.0, 9.0, 0.0 });
                           }));
        assert(bvp_data.mesh().nodes().group("fixed_nodes").node_indices.size() == 3);
        bvp_data.add_bc(DirichletBC("fixed_nodes", { 0, 0, BCApplicationType::NODE }));
        bvp_data.add_bc(DirichletBC("fixed_nodes", { 1, 0, BCApplicationType::NODE }));
        bvp_data.add_bc(DirichletBC("fixed_nodes", { 2, 0, BCApplicationType::NODE }));

        bvp.solve_time_steps(bvp_data);
    }

    return 0;
}