#include "../UnitTestHelpers.h"

#include "fea_analysis/FractureAnalysis.h"
#include "flex_mesh/algorithms/Creation.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/UnitElements.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "fea_io/common/MPIUtils.h"
#include "sparse_linear_algebra/LinAlgFactory.h"

#include "constitutive/enrichment/properties/TractionFreeInterfaceProperties.h"

using namespace fea_analysis;

int main(int argc, char ** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::init(
        sparse_linear_algebra::LinAlgFactory::LinAlgBackend::MKL, &argc, &argv);

    // A vector of { 0, 1, ..., 10 }
    std::vector<double> p(11);
    std::iota(p.begin(), p.end(), 0);
    // A vector of { 0, 1, ..., 15 }
    std::vector<double> z(16);
    std::iota(z.begin(), z.end(), 0);

    // Creates grid
    auto mesh_ptr = std::make_unique<Mesh>();
    *mesh_ptr = flex_mesh::make_grid<Mesh>(p, p, z, false);
    flex_mesh::write_vtk(*mesh_ptr, (get_output_directory(__FILE__) / "mesh.vtk").string());

	// Create a constitutive relation
    const double E = 1e9;
    const double nu = 0.3;
    constitutive::LinearElasticProperties properties;
    properties.set_isotropic_E_nu(E, nu);
    properties.set_thermal_expansion_isotropic(1e-4);

    FractureAnalysis analysis;
    auto weak_mesh_ptr = mesh_ptr.get();
    analysis.set_input(
        std::move(mesh_ptr),
        get_test_name(__FILE__),
        (get_output_directory(__FILE__) / "output").string());
    analysis.set_field_model(
        finite_element::BasalFieldSelection(),
        "InfinitesimalStrain",
        properties);
    for (auto& bc : get_unit_tension_z_bcs(*weak_mesh_ptr))
    {
        analysis.add_boundary_condition(std::move(bc));
    }
    analysis.setup_analysis();

    // Insert a spherical enrichment surface defined by the radius and center below.
    constitutive::TractionFreeInterfaceProperties enrichment_properties;
    const double sphere_radius1 = 4.0;
    const simple_math::Point sphere_center1 = { 5.0, 5.0, 5.0 };
    const double sphere_radius2 = 4.0;
    const simple_math::Point sphere_center2 = { 5.0, 5.0, 10.0 };
    analysis.insert_precracks_by_function_with_boolean(
        {
            {
                [sphere_radius = sphere_radius1, sphere_center = sphere_center1](const simple_math::Point& p) -> const double
                {
                    const double signed_dist = simple_math::Vector3D(
                        p[0] - sphere_center[0],
                        p[1] - sphere_center[1],
                        p[2] - sphere_center[2]).norm() - sphere_radius;
                    return signed_dist;
                },
                enrichment_properties,
                nullptr
            },
            {
                [sphere_radius = sphere_radius2, sphere_center = sphere_center2](const simple_math::Point& p) -> const double
                {
                    const double signed_dist = simple_math::Vector3D(
                        p[0] - sphere_center[0],
                        p[1] - sphere_center[1],
                        p[2] - sphere_center[2]).norm() - sphere_radius;
                    return signed_dist;
                },
                enrichment_properties,
                nullptr
            }
        });

    analysis.solve_time_steps();

	return 0;
}