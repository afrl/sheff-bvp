#include "../UnitTestHelpers.h"
#include "bvp/bvp/LinearElasticBVP.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "mesh/algorithms/Creation.h"
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
    bvp.set_input(mesh::make_mixed_element_grid({ 0, 0, 0 }, { 1, 1, 1 }, { 3, 3, 3 }),
                  get_test_name(__FILE__),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);
    bvp.set_field_models(finite_element::BasalFieldSelection(),
                         "InfinitesimalStrain",
                         "isotropic",
                         bvp_data);
    add_unit_tension_x_bcs(bvp_data);
    bvp.solve_time_steps(bvp_data);

    std::vector<finite_element::QuantityCoordSysPair<finite_element::VolumeQuantity>> flags = {
        { finite_element::VolumeQuantity::Strain, finite_element::CoordinateSystem::Local },
        { finite_element::VolumeQuantity::Stress, finite_element::CoordinateSystem::Local }
    };
    finite_element::MultiQuantityCalculationFlags<finite_element::VolumeQuantity> calc_flags(flags,
                                                                                             flags,
                                                                                             flags);

    auto element_workspace_ptr = std::make_unique<finite_element::ElementWorkspace>();
    simple_math::GeneralVector elem_solution_vector;
    for (auto& fe_accessor : bvp_data.fe_mesh().finite_elements())
    {
        if (mesh::get_category(fe_accessor.element.type()) == mesh::ElementTypeCategory::SOLID)
        {
            // Provide the current solution in case the element requries it.
            elem_solution_vector =
              assembler::extract_element_solution_vector(bvp_data.solution(),
                                                         bvp_data.dof_manager().get_dof_maps(),
                                                         fe_accessor.element.local_index());

            auto params = finite_element::FiniteElement::make_params(fe_accessor.element,
                                                                     elem_solution_vector,
                                                                     fe_accessor.layout,
                                                                     *element_workspace_ptr);

            auto&& result = finite_element::FiniteElement::calculate_output(params, calc_flags);

            auto& quad_local_strain = result.quad_results[0];
            auto& quad_local_stress = result.quad_results[1];
            for (auto q : utility::range(quad_local_strain.rows()))
            {
                assert(fabs(quad_local_strain(q, 0) - 1.0) < 1e-8);
                assert(fabs(quad_local_strain(q, 1) + 0.3) < 1e-8);
                assert(fabs(quad_local_strain(q, 2) + 0.3) < 1e-8);
                assert(fabs(quad_local_strain(q, 3)) < 1e-8);
                assert(fabs(quad_local_strain(q, 4)) < 1e-8);
                assert(fabs(quad_local_strain(q, 5)) < 1e-8);

                assert(fabs(quad_local_stress(q, 0) - E) < 1e1);
                for (auto i : utility::range(1, 6))
                {
                    assert(fabs(quad_local_stress(q, i)) < 1e1);
                }
            }

            auto& node_local_strain = result.node_results[0];
            auto& node_local_stress = result.node_results[1];
            for (auto i : utility::range(node_local_strain.rows()))
            {
                assert(fabs(node_local_strain(i, 0) - 1.0) < 1e-8);
                assert(fabs(node_local_strain(i, 1) + 0.3) < 1e-8);
                assert(fabs(node_local_strain(i, 2) + 0.3) < 1e-8);
                assert(fabs(node_local_strain(i, 3)) < 1e-8);
                assert(fabs(node_local_strain(i, 4)) < 1e-8);
                assert(fabs(node_local_strain(i, 5)) < 1e-8);

                assert(fabs(node_local_stress(i, 0) - E) < 1e1);
                for (auto i : utility::range(1, 6))
                {
                    assert(fabs(node_local_stress(i, i)) < 1e1);
                }
            }
        }
    }

    return 0;
}
