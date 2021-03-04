#include "../UnitTestHelpers.h"
#include "bvp/bvp/FractureBVP.h"
#include "constitutive/enrichment/deformation/bonded/BondedInterfaceProperties.h"
#include "constitutive/enrichment/deformation/traction_free/TractionFreeInterfaceProperties.h"
#include "constitutive/enrichment/direction_criteria/PositionFunctionDirectionCriterion.h"
#include "constitutive/enrichment/growth_criteria/AlwaysGrowGrowthCriterion.h"
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

    bvp::BVPGlobalOptions::get().save_global_matrices = true;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    constitutive::add_field_properties(
      "isotropic",
      constitutive::LinearElasticProperties::create_isotropic_E_nu(1e9, 0.3));

    auto surf_props = constitutive::add_enrichment_properties(
      "bonded",
      constitutive::BondedInterfaceProperties::create(1e20));

    auto crack_0_dir_crit = constitutive::add_direction_criterion(
      "crack_0_dir_crit",
      constitutive::PositionFunctionDirectionCriterion(
        [](const int& crack_index, const simple_math::Point& p) -> simple_math::Vector3D {
            return simple_math::Vector3D(0, 0, 1);
        }));

    auto growth_criterion =
      constitutive::ensure_growth_criterion("always_grow",
                                            constitutive::AlwaysGrowGrowthCriterion());

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Define boundary value problem

    auto bvp = FractureBVP::get();
    BVPData bvp_data;
    bvp.set_input(mesh::make_uniform_linear_grid({ 2, 2, 2 }, { 3, 3, 3 }, { 1, 1, 1 }),
                  get_test_name(__FILE__),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);
    bvp.set_field_models(finite_element::BasalFieldSelection(),
                         "InfinitesimalStrain",
                         "isotropic",
                         bvp_data);
    // add_unit_tension_z_bcs(bvp_data);

    // 1) Insert a spherical enrichment surface defined by the radius and center below.
    // 2) Insert a curved enrichment surface defined by the radius and center below that
    //    is initiated outside the sphere.
    /*
    const double sphere_radius             = 3.7;
    const simple_math::Point sphere_center = { 0, 0, 0 };
    constitutive::BondedInterfaceProperties enrichment_properties;
    bvp.insert_precracks_by_function_outside_previous(
      bvp_data,
      { { // surface 0
          [sphere_radius, sphere_center](const simple_math::Point& p) -> const double {
              const double signed_dist = simple_math::Vector3D(p[0] - sphere_center[0],
                                                               p[1] - sphere_center[1],
                                                               p[2] - sphere_center[2])
                                           .norm() -
                                         sphere_radius;
              return signed_dist;
          },
          surf_props,
          constitutive::FieldPropertiesCollectionRef() } });
          */


    constitutive::BondedInterfaceProperties enrichment_properties;
    bvp.insert_precracks(bvp_data,
                         { { // surface 0
                             { 2.5, 2.5, 2.5 },
                             surf_props,
                             growth_criterion,
                             crack_0_dir_crit } });

    bvp.solve_time_steps(bvp_data);

    auto& lin_alg = static_cast<LinearElasticLinearAlgebra&>(*bvp_data.gl_lin_alg);
    auto k        = lin_alg.stiffness->get_as_eigen_matrix();

    Eigen::EigenSolver<simple_math::GeneralMatrix> eigensolver;
    eigensolver.compute(k);
    auto eigenvalues  = eigensolver.eigenvalues();
    auto eigenvectors = eigensolver.eigenvectors();

    using EigenPair = std::tuple<double, simple_math::GeneralVector>;
    std::vector<EigenPair> eigs;
    for (auto i : utility::range(eigenvalues.size()))
    {
        eigs.push_back({ eigenvalues(i).real(), eigenvectors.col(i).real() });
    }
    // Sort eigs into descending order based on real part of eigenvalues
    std::sort(eigs.begin(), eigs.end(), [](const EigenPair& lhs, const EigenPair& rhs) -> bool {
        return std::get<0>(lhs) > std::get<0>(rhs);
    });

    mesh::Mesh m;
    auto n_basal_fields = bvp_data.fe_mesh().get_finite_element(0).layout.get_n_basal_fields();
    for (auto i : utility::range(n_basal_fields))
    {
        m.copy(bvp_data.mesh());
        m.elements().back().global_index() = i;
    }

    for (auto mode : utility::range(eigenvalues.size()))
    {
        auto& eig = eigs[mode];
        std::cout << "Writing mode " << mode << " with eigenvalue " << std::get<0>(eig)
                  << std::endl;
        Eigen::Map<simple_math::GeneralMatrix> u(std::get<1>(eig).data(), m.nodes().size(), 3);
        // Write the eigenvalue as cell data so you can visualize it
        std::vector<double> cell_eigenvalue(m.elements().size(), std::get<0>(eig));
        mesh::write_vtk(m,
                        "mode_" + std::to_string(mode) + ".vtk",
                        false,
                        {},
                        cell_eigenvalue,
                        {},
                        {},
                        u);
    }

    return 0;
}