#include "fea_analysis/LinearAnalysis.h"

#include "flex_mesh/algorithms/UnitElements.h"
#include "constitutive/elasticity/LinearElasticityConstitutiveRelation.h"
#include "fea_io/common/MPIUtils.h"

using namespace fea_analysis;

int main(int argc, char ** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);

	// Use a single linear hex for now
	auto mesh_ptr = flex_mesh::get_unit_linear_hex_ptr<fea_analysis::Mesh>();

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
        [](const fea_analysis::Node& n) -> bool {
        return n.x() == -1.0;
    });
	auto fix_u_at_x_neg_1 = std::unique_ptr<DirchletBC>(new DirchletBC(x_neg_1_selection, { 0 , 0 }));
    // Fix v along the y = -1 side
    NodeSelection y_neg_1_selection;
    y_neg_1_selection.add_criteria(
        [](const fea_analysis::Node& n) -> bool {
        return n.y() == -1.0;
    });
    auto fix_v_at_y_neg_1 = std::unique_ptr<DirchletBC>(new DirchletBC(y_neg_1_selection, { 1 , 0 }));
    // Fix w along the z = -1 side
    NodeSelection z_neg_1_selection;
    z_neg_1_selection.add_criteria(
        [](const fea_analysis::Node& n) -> bool {
        return n.z() == -1.0;
    });
    auto fix_w_at_z_neg_1 = std::unique_ptr<DirchletBC>(new DirchletBC(z_neg_1_selection, { 2 , 0 }));
	// Specify u along the x = 1 side
    NodeSelection x_pos_1_selection;
    x_pos_1_selection.add_criteria(
        [](const fea_analysis::Node& n) -> bool {
        return n.x() == 1.0;
    });
	auto fix_u_at_x_pos_1 = std::unique_ptr<DirchletBC>(new DirchletBC(x_pos_1_selection, { 0 , 1 }));

	LinearElasticAnalysis analysis;
	analysis.set_mesh(std::move(mesh_ptr));
	analysis.add_constitutive_relation("IsotropicMaterial", std::move(constitutive_relation_ptr));
	analysis.set_elements_constitutive_relation(fea_analysis::ElementSelection(), "IsotropicMaterial");
	analysis.add_boundary_condition(std::move(fix_u_at_x_neg_1));
    analysis.add_boundary_condition(std::move(fix_v_at_y_neg_1));
    analysis.add_boundary_condition(std::move(fix_w_at_z_neg_1));
	analysis.add_boundary_condition(std::move(fix_u_at_x_pos_1));
    analysis.setup_analysis();
	analysis.solve_time_steps();

	return 0;
}