#include "../UnitTestHelpers.h"

#include "fea_analysis/FractureAnalysis.h"
#include "flex_mesh/algorithms/Creation.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/UnitElements.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "fea_io/common/MPIUtils.h"
#include "sparse_linear_algebra/LinAlgFactory.h"

#include "constitutive/enrichment/properties/BondedInterfaceProperties.h"

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

    // Creates grid
    auto mesh_ptr = std::make_unique<Mesh>();
    *mesh_ptr = flex_mesh::make_grid<Mesh>(p, p, p, false);
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

    // Apply unit tension along z-axis
    for (auto& bc : get_unit_tension_z_bcs(*weak_mesh_ptr))
    {
        analysis.add_boundary_condition(std::move(bc));
    }

    analysis.setup_analysis();

    constitutive::LinearElasticProperties sphere_solid_properties;
    sphere_solid_properties.set_isotropic_E_nu(10e9, 0.3);
    sphere_solid_properties.set_thermal_expansion_isotropic(1e-4);

    // 1) Insert a spherical enrichment surface defined by the radius and center below.
    // 2) Insert a curved enrichment surface defined by the radius and center below that
    //    is initiated outside the sphere.
    const double sphere_radius = 3.7;
    const simple_math::Point sphere_center = { 5.0, 5.0, 5.0 };
    constitutive::BondedInterfaceProperties enrichment_properties;
    std::vector<SurfaceInsertionByLevelSet> insertions = { 
        { // surface 0
            [sphere_radius, sphere_center](const simple_math::Point& p) -> const double
            {
                const double signed_dist = simple_math::Vector3D(
                    p[0] - sphere_center[0],
                    p[1] - sphere_center[1],
                    p[2] - sphere_center[2]).norm() - sphere_radius;
                return signed_dist;
            },
            enrichment_properties,
            &sphere_solid_properties
        }};
    analysis.insert_precracks_by_function(insertions);

    analysis.solve_time_steps();

	return 0;
}