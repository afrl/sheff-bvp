#include "../UnitTestHelpers.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "bvp/bvp/LinearElasticBVP.h"
#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Partitioning.h"
#include "mesh/algorithms/UnitElements.h"
#include <chrono>
#include <thread>

using namespace bvp;

int main(int argc, char ** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::get().init(
      sparse_linear_algebra::LinAlgFactory::LinAlgBackend::PETSC,
      {});

    auto n_ranks = utility::MPIUtils::get_num_local_ranks();
    utility::MPIUtils::print_process_ids();
    std::this_thread::sleep_for(std::chrono::seconds(20));

    // A 3x3x1 grid
    auto total_mesh = mesh::make_linear_grid({ 0, 1. / 3., 2. / 3., 1.0 },
                                             { 0, 1. / 3., 2. / 3., 1.0 },
                                             { 0, 1.0 });
    // Delete the middle element to leave a hole
    total_mesh.elements().remove_stable(0);

    mesh::partition_metis(total_mesh, utility::MPIUtils::get_num_local_ranks());
    if (utility::MPIUtils::is_root())
    {
        mesh::write_fea_io(total_mesh, "./", "partitioned_mesh", "partitioned_mesh");
    }

    auto mesh_part = mesh::extract_partition_mesh(total_mesh, utility::MPIUtils::get_local_rank());
    mesh::write_vtk(mesh_part,
                    "mesh_part_" + std::to_string(utility::MPIUtils::get_local_rank()) + ".vtk");

    const double E  = 1e9;
    const double nu = 0.3;
    constitutive::add_field_properties(
      "isotropic",
      constitutive::LinearElasticProperties::create_isotropic_E_nu(E, nu));

    auto bvp = LinearElasticBVP::get();
    BVPData bvp_data;
    bvp.set_input(std::move(mesh_part),
                  get_test_name(__FILE__),
                  (get_output_directory(__FILE__) / "output").string(),
                  bvp_data);
    bvp.set_field_models(finite_element::BasalFieldSelection(),
                         "InfinitesimalStrain",
                         "isotropic",
                         bvp_data);
    add_unit_tension_x_bcs(bvp_data);
    bvp.solve_time_steps(bvp_data);

    return 0;
}
