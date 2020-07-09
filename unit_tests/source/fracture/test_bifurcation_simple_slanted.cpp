#include "../UnitTestHelpers.h"

#include "fea_analysis/FractureAnalysis.h"
#include "flex_mesh/algorithms/Creation.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/UnitElements.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "constitutive/enrichment/properties/TractionFreeInterfaceProperties.h"
#include "constitutive/enrichment/direction_criteria/PositionFunctionDirectionCriterion.h"
#include "constitutive/enrichment/growth_criteria/PositionFunctionGrowthCriterion.h"
#include "fea_io/common/MPIUtils.h"
#include "sparse_linear_algebra/LinAlgFactory.h"

using namespace fea_analysis;

int main(int argc, char ** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::init(
        sparse_linear_algebra::LinAlgFactory::LinAlgBackend::AUTO, &argc, &argv);

    // Creates a NxNxN grid from 0 to 10
    auto mesh_ptr = std::make_unique<Mesh>();
	*mesh_ptr = flex_mesh::make_grid<Mesh>({ 0, 1, 2, 3 }, { 0, 1, 2 }, { 0, 1, 2, 3 }, false);
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
    for (auto& bc : get_unit_tension_x_z_bcs(*weak_mesh_ptr))
    {
        analysis.add_boundary_condition(std::move(bc));
    }

    int phase = 0;
    constitutive::PositionFunctionDirectionCriterion direction_criterion(
        [](const int& crack_index, const simple_math::Point& p) -> simple_math::Vector3D
    {
		if (crack_index == 0)
		{
			if (p[1] <= 1)
			{
				return simple_math::Vector3D(0, 0, 1);
			}
			else if (p[0] > 1.5)
			{
				return simple_math::Vector3D(0, 0.3, 1);
			}
			else
			{
				return simple_math::Vector3D(0, -0.3, 1);
			}
		}
		else
		{
			if (p[1] >= 1)
			{
				return simple_math::Vector3D(1, 0, 0);
			}
			else if (p[2] < 1.5)
			{
				return simple_math::Vector3D(1, 0.3, 0);
			}
			else
			{
				return simple_math::Vector3D(1, -0.3, 0);
			}
		}
    });
	constitutive::PositionFunctionGrowthCriterion growth_criterion(
		[&phase](const int& crack_index, const simple_math::Point& p) -> bool
	{
        switch (phase)
        {
        case 0:
            return false;
        case 1:
            return crack_index == 0;
        case 2:
            return crack_index == 1;
        }
	});
    constitutive::TractionFreeInterfaceProperties enrichment_properties;

    analysis.setup_analysis();
    analysis.insert_precracks({
        {
            { 1.5, 0.2, 1.5 }, // point
            enrichment_properties,
            growth_criterion,
            direction_criterion
        },
        {
            { 1.5, 1.8, 1.5 }, // point
            enrichment_properties,
            growth_criterion,
            direction_criterion
        }
    });

    // Let the cracks grow past each other
    phase = 1;
    analysis.try_to_grow_precracks();
    phase = 2;
    analysis.try_to_grow_precracks();

    analysis.solve_time_steps();

	return 0;
}