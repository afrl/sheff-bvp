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

	auto N = 9;
    std::vector<double> grid_offsets(N);
    for (auto i : utility::indices(grid_offsets))
    {
        grid_offsets[i] = 2.0 * double(i) / double(grid_offsets.size() - 1) - 1.0;
    }

    // Creates a NxNxN grid from -1 to 1
    auto mesh_ptr = std::make_unique<Mesh>();
    *mesh_ptr = flex_mesh::make_grid<Mesh>(
        grid_offsets, grid_offsets, { 0, 0.2, 0.4 }, false);
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
    for (auto& bc : get_unit_tension_x_bcs())
    {
        analysis.add_boundary_condition(std::move(bc));
    }

    constitutive::PositionFunctionDirectionCriterion direction_criterion(
        [](const int& crack_index, const simple_math::Point& point) -> simple_math::Vector3D
    {
        return simple_math::Vector3D(1, 0, 0);
    });
    for (auto& e_ptr : weak_mesh_ptr->elements())
    {
        e_ptr->set_direction_criterion(direction_criterion);
    }

    analysis.setup_analysis();
    simple_math::Point initial_coord = { 0.525, 0.2, 0.2 };
    //analysis.insert_precracks({ initial_coord });
    analysis.solve_time_steps();

	return 0;
}