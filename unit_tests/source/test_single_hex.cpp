#include "UnitTestHelpers.h"

#include "fea_analysis/LinearAnalysis.h"
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

    // Use a single linear hex
    auto mesh_ptr = flex_mesh::get_unit_linear_hex_ptr<fea_analysis::Mesh>();

    // Create a constitutive relation
    const double E = 1e9;
    const double nu = 0.3;
    auto properties = std::make_unique<constitutive::LinearElasticProperties>();
    properties->set_isotropic_E_nu(E, nu);
    properties->set_thermal_expansion_isotropic(1e-4);
    auto constitutive_relation_ptr = std::make_unique<constitutive::LinearElasticityConstitutiveRelation>();
    constitutive_relation_ptr->set_elastic_properties(std::move(properties));

    LinearElasticAnalysis analysis;
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
    analysis.setup_analysis();
    analysis.solve_time_steps();

    return 0;
}
