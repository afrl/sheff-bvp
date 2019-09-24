#include "fea_analysis/FractureAnalysis.h"

#include "flex_mesh/algorithms/Creation.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/UnitElements.h"

#include "finite_element/elasticity/InfinitesimalStrain.h"

#include "dof_manager/DofMap.h"

#include "constitutive/elasticity/LinearElasticityConstitutiveRelation.h"

#include "fea_io/common/MPIUtils.h"

#include "simple_math/Utility.h"

using namespace fea_analysis;

int main(int argc, char ** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);

    auto N = 9;
    std::vector<double> grid_offsets(N);
    for (auto i : utility::indices(grid_offsets))
    {
        grid_offsets[i] = double(i) / double(grid_offsets.size() - 1);
    }

    // Creates a NxNxN grid from 0 to 1
    auto mesh_ptr = std::make_unique<Mesh>();
    *mesh_ptr = flex_mesh::make_grid<Mesh>(
        grid_offsets, grid_offsets, { 0, 0.1 }, false);

    //flex_mesh::write_vtk(*mesh_ptr, "./SelfIntersectGrid.vtk", false, false);

    auto properties = std::make_unique<constitutive::LinearElasticProperties>();
    properties->set_isotropic_E_nu(1e9, 0.3);
    properties->set_thermal_expansion_isotropic(1e-4);
    auto constitutive_relation_ptr = std::make_unique<constitutive::LinearElasticityConstitutiveRelation>();
    constitutive_relation_ptr->set_elastic_properties(std::move(properties));

	// Create the boundary conditions
	// Fix u along the x = -1 side
    NodeSelection x_0_selection;
    x_0_selection.add_criteria(
        [](const fea_analysis::Node& n) -> bool {
        return n.x() == 0.0;
    });
	auto fix_u_at_x_0 = std::unique_ptr<DirchletBC>(new DirchletBC(x_0_selection, { 0 , 0 }));
    // Fix v along the y = -1 side
    NodeSelection y_0_selection;
    y_0_selection.add_criteria(
        [](const fea_analysis::Node& n) -> bool {
        return n.y() == 0.0;
    });
    auto fix_v_at_y_0 = std::unique_ptr<DirchletBC>(new DirchletBC(y_0_selection, { 1 , 0 }));
    // Fix w along the z = -1 side
    NodeSelection z_0_selection;
    z_0_selection.add_criteria(
        [](const fea_analysis::Node& n) -> bool {
        return n.z() == 0.0;
    });
    auto fix_w_at_z_0 = std::unique_ptr<DirchletBC>(new DirchletBC(z_0_selection, { 2 , 0 }));
	// Specify v along the y = 1 side
    NodeSelection y_pos_1_selection;
    y_pos_1_selection.add_criteria(
        [](const fea_analysis::Node& n) -> bool {
        return n.y() == 1.0;
    });
	auto fix_v_at_y_pos_1 = std::unique_ptr<DirchletBC>(new DirchletBC(y_pos_1_selection, { 1 , 1 }));

    FractureAnalysis analysis;
    auto weak_mesh_ptr = mesh_ptr.get();
	analysis.set_mesh(std::move(mesh_ptr));
	analysis.add_constitutive_relation("IsotropicMaterial", std::move(constitutive_relation_ptr));
	analysis.set_elements_constitutive_relation(fea_analysis::ElementSelection(), "IsotropicMaterial");
	analysis.add_boundary_condition(std::move(fix_u_at_x_0));
    analysis.add_boundary_condition(std::move(fix_v_at_y_0));
    analysis.add_boundary_condition(std::move(fix_w_at_z_0));
	analysis.add_boundary_condition(std::move(fix_v_at_y_pos_1));
	
    bool entered_circle_region = false;
    constitutive::PositionFunctionDirectionCriterion direction_criterion(
        [&entered_circle_region](const int & crack_index, const simple_math::Point & point) -> simple_math::Vector3D
    {
        simple_math::Vector3D crack_normal;
        if (point[0] <= 0.5 && !entered_circle_region)
        {
            crack_normal << 0, -1, 0;
        }
        else if (point[0] <= 0.3 && point[1] <= 0.5 && entered_circle_region)
        {
            crack_normal << -1, 0, 0;
        }
        else
        {
            auto point_vec = simple_math::point_as_vector(point);
            simple_math::Vector3D point_2D(point_vec[0] - 0.5, point_vec[1] - 0.5, 0);
            crack_normal = point_2D.normalized();
            entered_circle_region = true;
        }
        return crack_normal;
    });
    for (auto& e_ptr : weak_mesh_ptr->elements())
    {
        e_ptr->set_direction_criterion(direction_criterion);
    }
    
    analysis.setup_analysis();
    simple_math::Point initial_coord = { 0.075, 0.2, 0 };
    analysis.insert_precracks({ initial_coord });
    analysis.solve_time_steps();

	return 0;
}