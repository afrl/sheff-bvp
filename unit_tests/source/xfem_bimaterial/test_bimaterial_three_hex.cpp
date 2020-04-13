#include "../UnitTestHelpers.h"

#include "fea_analysis/FractureAnalysis.h"
#include "flex_mesh/algorithms/Creation.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/UnitElements.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "fea_io/common/MPIUtils.h"
#include "sparse_linear_algebra/LinAlgFactory.h"


#include "constitutive/enrichment/direction_criteria/PositionFunctionDirectionCriterion.h"
#include "constitutive/enrichment/growth_criteria/AlwaysGrowGrowthCriterion.h"
#include "constitutive/enrichment/properties/BondedInterfaceProperties.h"

using namespace fea_analysis;

int main(int argc, char ** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::init(
        sparse_linear_algebra::LinAlgFactory::LinAlgBackend::MKL, &argc, &argv);

    // Creates grid
    auto mesh_ptr = std::make_unique<Mesh>();
    *mesh_ptr = flex_mesh::make_grid<Mesh>(
        { 0, 1, 2, 3 }, { 0, 1 }, { 0, 1 }, false);
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
    for (auto& bc : get_unit_tension_x_bcs(*weak_mesh_ptr))
    {
        analysis.add_boundary_condition(std::move(bc));
    }

    constitutive::AlwaysGrowGrowthCriterion growth_criterion;
    constitutive::PositionFunctionDirectionCriterion direction_criterion(
        [](const int& crack_index, const simple_math::Point& point) -> simple_math::Vector3D
    {
        return simple_math::Vector3D(1, 0, 0);
    });
    constitutive::BondedInterfaceProperties enrichment_properties;

    analysis.setup_analysis();
    simple_math::Point initial_coord = { 1.525, 0.2, 0.2 };
    analysis.insert_precracks({
        { // surface 0
                initial_coord, // initation_point
                enrichment_properties,
                growth_criterion,
                direction_criterion
        } });
    analysis.solve_time_steps();

	return 0;
}