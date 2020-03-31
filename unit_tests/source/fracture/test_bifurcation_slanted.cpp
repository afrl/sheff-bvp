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

    // Creates a NxNxN grid from 0 to 10
    auto mesh_ptr = std::make_unique<Mesh>();
	*mesh_ptr = flex_mesh::make_grid<Mesh>(p, p, p, false);
    flex_mesh::write_vtk(*mesh_ptr, (get_output_directory(__FILE__) / "mesh.vtk").string());

	// Create a constitutive relation
    const double E = 1e9;
    const double nu = 0.3;
    auto properties = std::make_unique<constitutive::LinearElasticProperties>();
    properties->set_isotropic_E_nu(E, nu);
    properties->set_thermal_expansion_isotropic(1e-4);
	auto constitutive_relation_ptr = std::make_unique<constitutive::LinearElasticityConstitutiveRelation>();
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
    for (auto& bc : get_unit_tension_x_z_bcs(*weak_mesh_ptr))
    {
        analysis.add_boundary_condition(std::move(bc));
    }

    constitutive::PositionFunctionDirectionCriterion direction_criterion(
        [](const int& crack_index, const simple_math::Point& p) -> simple_math::Vector3D
    {
		if (crack_index == 0)
		{
			if (p[1] <= 5)
			{
				return simple_math::Vector3D(0, 0, 1);
			}
			else if (p[0] > 4.5)
			{
				return simple_math::Vector3D(0, 0.6, 1);
			}
			else
			{
				return simple_math::Vector3D(0, -0.6, 1);
			}
		}
		else
		{
			if (p[1] >= 5)
			{
				return simple_math::Vector3D(1, 0, 0);
			}
            // Note: for the two cases below, 0.6 instead of 0.45 leads to a
            // problem where one node is tied to the other side of the surface.
            // This needs investigation.
			else if (p[2] < 4.5)
			{
				return simple_math::Vector3D(1, 0.45, 0);
			}
			else
			{
				return simple_math::Vector3D(1, -0.45, 0);
			}
		}
    });
	constitutive::PositionFunctionPropagationCriterion propagation_criterion(
		[](const int& crack_index, const simple_math::Point& p) -> bool
	{
		return (p[1] < 4 && crack_index == 0)
			|| (p[1] > 6 && crack_index == 1);
	});

    for (auto& e_ptr : weak_mesh_ptr->elements())
    {
        e_ptr->set_direction_criterion(direction_criterion);
		e_ptr->set_propagation_criterion(propagation_criterion);
    }

    analysis.setup_analysis();
    analysis.insert_precracks({ 
		{ 4.5, 0.2, 4.5 },
		{ 4.5, 9.8, 4.5 }
	});

    // Let the cracks grow past each other
    constitutive::PositionFunctionPropagationCriterion propagation_criterion2(
        [](const int& crack_index, const simple_math::Point& p) -> bool
    {
        return crack_index == 0;
    });
    for (auto& e_ptr : weak_mesh_ptr->elements())
    {
        e_ptr->set_propagation_criterion(propagation_criterion2);
    }
    analysis.try_to_grow_precracks();
    
    constitutive::PositionFunctionPropagationCriterion propagation_criterion3(
        [](const int& crack_index, const simple_math::Point& p) -> bool
    {
        return crack_index == 1;
    });
    for (auto& e_ptr : weak_mesh_ptr->elements())
    {
        e_ptr->set_propagation_criterion(propagation_criterion3);
    }
    analysis.try_to_grow_precracks();

    analysis.solve_time_steps();

	return 0;
}