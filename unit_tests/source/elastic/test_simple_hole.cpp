#include "../UnitTestHelpers.h"

#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"

#include "fea_analysis/LinearAnalysis.h"

#include "fea_io/common/MPIUtils.h"

#include "flex_mesh/algorithms/Creation.h"
#include "flex_mesh/algorithms/Partitioning.h"

#include "sparse_linear_algebra/LinAlgFactory.h"

using namespace fea_analysis;

int main(int argc, char ** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::init(
        sparse_linear_algebra::LinAlgFactory::LinAlgBackend::MKL, &argc, &argv);

    utility::MPIUtils::print_process_ids();

    // Use a single linear hex
    auto total_mesh = flex_mesh::make_grid<fea_analysis::Mesh>(
            { 0, 1. / 3., 2. / 3., 1.0 },
            { 0, 1. / 3., 2. / 3., 1.0 },
            { 0, 1.0 },
            false);
    // Delete the middle element to leave a hole
    total_mesh.elements().remove(total_mesh.elements().begin() + 4);
    total_mesh.update_global_numbering();

    flex_mesh::partition_metis(total_mesh, utility::MPIUtils::get_num_local_ranks());

    auto mesh_ptr = std::make_unique<fea_analysis::Mesh>(
        flex_mesh::extract_partition_mesh(total_mesh, utility::MPIUtils::get_local_rank()));
    auto weak_mesh_ptr = mesh_ptr.get();

    // Create a constitutive relation
    const double E = 1e9;
    const double nu = 0.3;
    constitutive::LinearElasticProperties properties;
    properties.set_isotropic_E_nu(E, nu);
    properties.set_thermal_expansion_isotropic(1e-4);

    LinearElasticAnalysis analysis;
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
    analysis.setup_analysis();
    analysis.solve_time_steps();

    return 0;
}
