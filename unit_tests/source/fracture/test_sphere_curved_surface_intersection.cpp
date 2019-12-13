#include "../UnitTestHelpers.h"

#include "fea_analysis/FractureAnalysis.h"
#include "flex_mesh/algorithms/Creation.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/UnitElements.h"
#include "constitutive/elasticity/LinearElasticityConstitutiveRelation.h"
#include "fea_io/common/MPIUtils.h"
#include "sparse_linear_algebra/LinAlgFactory.h"

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
    auto properties = std::make_unique<constitutive::LinearElasticProperties>();
    properties->set_isotropic_E_nu(E, nu);
    properties->set_thermal_expansion_isotropic(1e-4);
	auto constitutive_relation_ptr =
        std::make_unique<constitutive::LinearElasticityConstitutiveRelation>();
	constitutive_relation_ptr->set_elastic_properties(std::move(properties));

    FractureAnalysis analysis;
    auto weak_mesh_ptr = mesh_ptr.get();
    analysis.set_input(
        std::move(mesh_ptr),
        get_test_name(__FILE__),
        (get_output_directory(__FILE__) / "output").string());
    analysis.add_constitutive_relation(
        "IsotropicMaterial",
        std::move(constitutive_relation_ptr));
    analysis.set_elements_constitutive_relation(
        fea_analysis::ElementSelection(),
        "IsotropicMaterial");

    // Apply unit tension along z-axis
    for (auto& bc : get_unit_tension_z_bcs(*weak_mesh_ptr))
    {
        analysis.add_boundary_condition(std::move(bc));
    }
    // Add an extra boundary condition to move the sphere in the center upwards
    {
        NodeSelection center_selection;
        center_selection.add_criteria(
            [](const fea_analysis::Node& n) -> bool {
            return n.x() == 5.0 && n.y() == 5.0 && n.z() == 5.0;
        });
        analysis.add_boundary_condition(std::unique_ptr<DirichletBC>(
            new DirichletBC(center_selection, { 0 , 0.0 })));
        analysis.add_boundary_condition(std::unique_ptr<DirichletBC>(
            new DirichletBC(center_selection, { 1 , 0.0 })));
        analysis.add_boundary_condition(std::unique_ptr<DirichletBC>(
            new DirichletBC(center_selection, { 2 , 0.5 })));
    }
    // Prevent rigid body rotation of the center
    {
        NodeSelection selection_1;
        selection_1.add_criteria(
            [](const fea_analysis::Node& n) -> bool {
            return n.x() == 6.0 && n.y() == 5.0 && n.z() == 5.0;
        });
        analysis.add_boundary_condition(std::unique_ptr<DirichletBC>(
            new DirichletBC(selection_1, { 0 , 0.0 })));
        analysis.add_boundary_condition(std::unique_ptr<DirichletBC>(
            new DirichletBC(selection_1, { 1 , 0.0 })));
        analysis.add_boundary_condition(std::unique_ptr<DirichletBC>(
            new DirichletBC(selection_1, { 2 , 0.5 })));
        NodeSelection selection_2;
        selection_2.add_criteria(
            [](const fea_analysis::Node& n) -> bool {
            return n.x() == 5.0 && n.y() == 6.0 && n.z() == 5.0;
        });
        analysis.add_boundary_condition(std::unique_ptr<DirichletBC>(
            new DirichletBC(selection_2, { 0 , 0.0 })));
        analysis.add_boundary_condition(std::unique_ptr<DirichletBC>(
            new DirichletBC(selection_2, { 1 , 0.0 })));
        analysis.add_boundary_condition(std::unique_ptr<DirichletBC>(
            new DirichletBC(selection_2, { 2 , 0.5 })));
    }

    analysis.setup_analysis();

    // 1) Insert a spherical enrichment surface defined by the radius and center below.
    // 2) Insert a curved enrichment surface defined by the radius and center below that
    //    is initiated outside the sphere.
    const double sphere_1_radius = 3.7;
    const simple_math::Point sphere_1_center = { 5.0, 5.0, 5.0 };
    const double sphere_2_radius = 8.0;
    simple_math::Point sphere_2_center = { 5, 5, 12 };
    analysis.insert_precracks_by_function_outside_previous(
        {
            [sphere_1_radius, sphere_1_center](const simple_math::Point& p) -> const double
            {
                const double signed_dist = simple_math::Vector3D(
                    p[0] - sphere_1_center[0],
                    p[1] - sphere_1_center[1],
                    p[2] - sphere_1_center[2]).norm() - sphere_1_radius;
                return signed_dist;
            },
            [sphere_2_radius, sphere_2_center](const simple_math::Point& p) -> const double
            {
                const double signed_dist = simple_math::Vector3D(
                    p[0] - sphere_2_center[0],
                    p[1] - sphere_2_center[1],
                    p[2] - sphere_2_center[2]).norm() - sphere_2_radius;
                return signed_dist;
            }
        });

    analysis.solve_time_steps();

	return 0;
}