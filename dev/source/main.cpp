#include "fea_analysis/LinearAnalysis.h"

#include "flex_mesh/algorithms/UnitElements.h"
#include "constitutive/elasticity/LinearElasticity.h"

using namespace fea_analysis;

int main()
{
	// Use a single linear hex for now
	auto mesh_ptr = flex_mesh::get_unit_linear_hex_ptr<fea_analysis::Mesh>();

	// Create a constitutive relation
	auto constitutive_relation_ptr = std::make_unique<constitutive::LinearElasticity>();
	const double E = 1e9;
	const double nu = 0.3;
	constitutive_relation_ptr->set_from_engineering_parameters(E, nu);

	// Create the boundary conditions
	// Fix u along the x = -1 side
	auto x_neg_1 = std::make_unique<flex_mesh::LambdaCriteria<fea_analysis::Node>>(
		[](const fea_analysis::Node& n) -> bool {
			return n.get_coord()[0] == -1.0;
		});
	auto x_neg_1_selection = std::make_unique<NodeSelection>();
	x_neg_1_selection->add_criteria(std::move(x_neg_1));
	auto fix_u_at_x_neg_1 = std::unique_ptr<DirchletBC>(new DirchletBC(std::move(x_neg_1_selection), { 0 , 0 }));
	// Specify u along the x = 1 side
	auto x_pos_1 = std::make_unique<flex_mesh::LambdaCriteria<fea_analysis::Node>>(
		[](const fea_analysis::Node& n) -> bool {
		return n.get_coord()[0] == 1.0;
	});
	auto x_pos_1_selection = std::make_unique<NodeSelection>();
	x_pos_1_selection->add_criteria(std::move(x_pos_1));
	auto fix_u_at_x_pos_1 = std::unique_ptr<DirchletBC>(new DirchletBC(std::move(x_pos_1_selection), { 0 , 1 }));

	LinearElasticAnalysis analysis;
	analysis.set_mesh(std::move(mesh_ptr));
	analysis.add_constitutive_relation("IsotropicMaterial", std::move(constitutive_relation_ptr));
	analysis.set_elements_constitutive_relation(fea_analysis::ElementSelection(), "IsotropicMaterial");
	analysis.add_boundary_condition(std::move(fix_u_at_x_neg_1));
	analysis.add_boundary_condition(std::move(fix_u_at_x_pos_1));
	analysis.solve_time_steps();

	return 0;
}