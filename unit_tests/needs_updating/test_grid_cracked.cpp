#include "bvp/FractureAnalysis.h"

#include "flex_mesh/algorithms/Creation.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/UnitElements.h"

#include "constitutive/elasticity/LinearElasticityConstitutiveRelation.h"
#include "fea_io/common/MPIUtils.h"

using namespace bvp;

int main(int argc, char ** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);

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

	// Create a constitutive relation
    const double E = 1e9;
    const double nu = 0.3;
    auto properties = std::make_unique<constitutive::LinearElasticProperties>();
    properties->set_isotropic_E_nu(E, nu);
    properties->set_thermal_expansion_isotropic(1e-4);
	auto constitutive_relation_ptr = std::make_unique<constitutive::LinearElasticityConstitutiveRelation>();
	constitutive_relation_ptr->set_elastic_properties(std::move(properties));

	// Create the boundary conditions
	// Fix u along the x = -1 side
    NodeSelection x_neg_1_selection;
	x_neg_1_selection.add_criteria(
        [](const bvp::Node& n) -> bool {
        return n.x() == -1.0;
    });
	auto fix_u_at_x_neg_1 = std::unique_ptr<DirchletBC>(new DirchletBC(x_neg_1_selection, { 0 , 0 }));
    // Fix v along the y = -1 side
    NodeSelection y_neg_1_selection;
    y_neg_1_selection.add_criteria(
        [](const bvp::Node& n) -> bool {
        return n.y() == -1.0;
    });
    auto fix_v_at_y_neg_1 = std::unique_ptr<DirchletBC>(new DirchletBC(y_neg_1_selection, { 1 , 0 }));
    // Fix w along the z = 0 side
    NodeSelection z_neg_selection;
    z_neg_selection.add_criteria(
        [](const bvp::Node& n) -> bool {
        return n.z() == 0.0;
    });
    auto fix_w_at_z_neg = std::unique_ptr<DirchletBC>(new DirchletBC(z_neg_selection, { 2 , 0 }));
	// Specify u along the x = 1 side
    NodeSelection x_pos_1_selection;
    x_pos_1_selection.add_criteria(
        [](const bvp::Node& n) -> bool {
        return n.x() == 1.0;
    });
	auto fix_u_at_x_pos_1 = std::unique_ptr<DirchletBC>(new DirchletBC(x_pos_1_selection, { 0 , 1 }));

    FractureAnalysis analysis;
    auto weak_mesh_ptr = mesh_ptr.get();
	analysis.set_mesh(std::move(mesh_ptr));
	analysis.add_constitutive_relation("IsotropicMaterial", std::move(constitutive_relation_ptr));
	analysis.set_elements_constitutive_relation(bvp::ElementSelection(), "IsotropicMaterial");
	analysis.add_boundary_condition(std::move(fix_u_at_x_neg_1));
    analysis.add_boundary_condition(std::move(fix_v_at_y_neg_1));
    analysis.add_boundary_condition(std::move(fix_w_at_z_neg));
	analysis.add_boundary_condition(std::move(fix_u_at_x_pos_1));

    bool entered_circle_region = false;
    constitutive::PositionFunctionDirectionCriterion direction_criterion(
        [&entered_circle_region](const int& crack_index, const simple_math::Point& point) -> simple_math::Vector3D
    {
        return simple_math::Vector3D(1, 0, 0);
    });
    for (auto& e_ptr : weak_mesh_ptr->elements())
    {
        e_ptr->set_direction_criterion(direction_criterion);
    }

    analysis.setup_analysis();
    simple_math::Point initial_coord = { 0.525, 0.2, 0.2 };
    analysis.insert_precracks({ initial_coord });
    analysis.solve_time_steps();

	return 0;
}