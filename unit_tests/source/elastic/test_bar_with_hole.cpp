#include "../UnitTestHelpers.h"

#include "fea_analysis/LinearAnalysis.h"
#include "flex_mesh/algorithms/UnitElements.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/Search.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "fea_io/common/MPIUtils.h"
#include "sparse_linear_algebra/LinAlgFactory.h"

using namespace fea_analysis;

int main(int argc, char ** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::init(
        sparse_linear_algebra::LinAlgFactory::LinAlgBackend::AUTO, &argc, &argv);
    init_output_directory(__FILE__);

    // Create a properties object
    const double E = 1e9;
    const double nu = 0.3;
    constitutive::LinearElasticProperties properties;
    properties.set_isotropic_E_nu(1e9, 0.3);
    properties.set_thermal_expansion_isotropic(1e-4);

    // Create the boundary conditions
    // Fix u along the x = -1 side
    NodeSelection x_neg_1_selection;
    x_neg_1_selection.add_criteria(
        [](const fea_analysis::Node& n) -> bool {
        return n.x() == 0.0;
    });
    auto fix_u_at_x_neg_1 = std::unique_ptr<DirichletBC>(new DirichletBC(x_neg_1_selection, { 0 , 0 }));
    auto fix_v_at_x_neg_1 = std::unique_ptr<DirichletBC>(new DirichletBC(x_neg_1_selection, { 1 , 0 }));
    auto fix_w_at_x_neg_1 = std::unique_ptr<DirichletBC>(new DirichletBC(x_neg_1_selection, { 2 , 0 }));
    // Specify u along the x = 1 side
    NodeSelection x_pos_1_selection;
    x_pos_1_selection.add_criteria(
        [](const fea_analysis::Node& n) -> bool {
        return n.x() == 100.0;
    });
    auto fix_u_at_x_pos_1 = std::unique_ptr<DirichletBC>(new DirichletBC(x_pos_1_selection, { 0 , 10 }));

    LinearElasticAnalysis analysis;
    analysis.set_input(
        "",
        get_test_data_filepath("bar_with_hole"),
        (get_output_directory(__FILE__) / "output").string());
    analysis.set_field_model(finite_element::BasalFieldSelection(), "InfinitesimalStrain", properties);
    analysis.add_boundary_condition(std::move(fix_u_at_x_neg_1));
    analysis.add_boundary_condition(std::move(fix_v_at_x_neg_1));
    analysis.add_boundary_condition(std::move(fix_w_at_x_neg_1));
    analysis.add_boundary_condition(std::move(fix_u_at_x_pos_1));
    analysis.setup_analysis();
    analysis.solve_time_steps();

    return 0;
}