#include "../UnitTestHelpers.h"
#include "bvp/bc/TractionBC.h"
#include "bvp/bvp/LinearElasticBVP.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "mesh/algorithms/UnitElements.h"

using namespace bvp;

int main(int argc, char** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::init(
      sparse_linear_algebra::LinAlgFactory::LinAlgBackend::AUTO,
      &argc,
      &argv);

    const double E  = 1e9;
    const double nu = 0.3;
    constitutive::add_field_properties(
      "isotropic",
      constitutive::LinearElasticProperties::create_isotropic_E_nu(E, nu));

    auto bvp = LinearElasticBVP::get();

    // Test traction for a linear hexahedron element
    {
        BVPData bvp_data;
        bvp.set_input(mesh::get_unit_element(mesh::ElementType::H3D8N),
                      get_test_name(__FILE__),
                      (get_output_directory(__FILE__) / "linear_output").string(),
                      bvp_data);
        bvp.set_field_models(finite_element::BasalFieldSelection(),
                             "InfinitesimalStrain",
                             "isotropic",
                             bvp_data);
        add_extent_node_groups(bvp_data.mesh());
        add_symmetry_about_0_planes_bcs(bvp_data);
        bvp_data.add_bc(TractionBC("pos_x", { { 1e6, 0, 0 } }));
        bvp.solve_time_steps(bvp_data);

        auto& constraint_matrix = *bvp_data.gl_lin_alg->constraint_matrices()[0];
        auto F                  = constraint_matrix.get_owned_nodal_neumann_terms();
        double force_sum        = 0.0;
        for (auto i : utility::range(F.size()))
        {
            force_sum += F[i];
            std::cout << i << "  " << F[i] << std::endl;
        }
        assert(fabs(force_sum - 4 * 1e6) < (4 * 1e6 * 1e-8));

        auto solution = bvp_data.solution().get_as_eigen_vector();
        for (auto i : bvp_data.mesh().nodes().group("pos_x").node_indices)
        {
            assert(fabs(solution(3 * i) - 2e-3) < 1e-8);
        }
    }

    // Test traction for a quadratic hexahedron element
    {
        BVPData bvp_data;
        bvp.set_input(mesh::get_unit_element(mesh::ElementType::H3D20N),
                      get_test_name(__FILE__),
                      (get_output_directory(__FILE__) / "quadratic_output").string(),
                      bvp_data);
        bvp.set_field_models(finite_element::BasalFieldSelection(),
                             "InfinitesimalStrain",
                             "isotropic",
                             bvp_data);
        add_extent_node_groups(bvp_data.mesh());
        add_symmetry_about_0_planes_bcs(bvp_data);
        bvp_data.add_bc(TractionBC("pos_x", { { 1e6, 0, 0 } }));
        bvp.solve_time_steps(bvp_data);

        auto& constraint_matrix = *bvp_data.gl_lin_alg->constraint_matrices()[0];
        auto F                  = constraint_matrix.get_owned_nodal_neumann_terms();
        double force_sum        = 0.0;
        for (auto i : utility::range(F.size()))
        {
            force_sum += F[i];
            std::cout << i << "  " << F[i] << std::endl;
        }
        assert(fabs(force_sum - 4 * 1e6) < (4 * 1e6 * 1e-8));

        auto solution = bvp_data.solution().get_as_eigen_vector();
        for (auto i : bvp_data.mesh().nodes().group("pos_x").node_indices)
        {
            assert(fabs(solution(3 * i) - 2e-3) < 1e-8);
        }
    }

    return 0;
}
