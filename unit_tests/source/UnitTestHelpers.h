#pragma once

#include "bvp/bc/DirchletBC.h"
#include "bvp/bvp/BVPData.h"
#include "fea_io/common/MPIUtils.h"
#include "mesh/algorithms/Groups.h"
#include "mesh/algorithms/Search.h"
#include "sparse_linear_algebra/LinAlgFactory.h"
#include "utility/Exceptions.h"
#include "utility/Filesystem.h"
#include "utility/MPIUtils.h"

namespace bvp
{
/// Given the file name for a test data file, this returns the full filepath
static const std::string get_test_data_filepath(const std::string& data_filename)
{
    // Assumes test data is one directory up in another directory called "data"
    auto this_file_path = fs::canonical(fs::path(__FILE__));
    auto data_path =
      this_file_path.parent_path().parent_path().append("data").append(data_filename);
    return data_path.string();
}

static std::string get_test_name(const std::string& test_filepath)
{
    return fs::path(test_filepath).stem().string();
}

static fs::path get_output_directory(const std::string& test_filepath)
{
    auto test_name = get_test_name(test_filepath);
    std::string output_dir;
    if (fs::exists(fs::path("./") / "Output"))
    {
        output_dir = "./Output/" + test_name;
    }
    else
    {
        output_dir = "./" + test_name;
    }
    if (!fs::exists(output_dir))
    {
        fs::create_directory(output_dir);
    }
    return fs::path(output_dir);
}


static void init_output_directory(const std::string& test_filepath)
{
    auto output_dir = get_output_directory(test_filepath);
    if (!fs::exists(output_dir))
    {
        fs::create_directory(output_dir);
    }
    else
    {
        if (!fs::is_directory(output_dir))
        {
            throw _UTIL_MSG_EXCEPT("The directory for output (" + output_dir.string() +
                                   ") already exists and is not a directory.");
        }
        else
        {
            fs::remove_all(output_dir);
            fs::create_directory(output_dir);
        }
    }
}

static void add_extent_node_groups(mesh::Mesh& mesh)
{
    if (mesh.nodes().group("neg_x").node_indices.size() > 0)
    {
        return;
    }

    auto extent     = mesh::find_global_extent(mesh);
    auto& min_coord = extent.first;
    auto& max_coord = extent.second;

    mesh::add_to_group(
      mesh,
      "neg_x",
      mesh::NodeSelection(
        [val = min_coord[0]](const mesh::NodeConstAccessor& n) -> bool { return n.x() == val; }));
    mesh::add_to_group(
      mesh,
      "pos_x",
      mesh::NodeSelection(
        [val = max_coord[0]](const mesh::NodeConstAccessor& n) -> bool { return n.x() == val; }));
    mesh::add_to_group(
      mesh,
      "neg_y",
      mesh::NodeSelection(
        [val = min_coord[1]](const mesh::NodeConstAccessor& n) -> bool { return n.y() == val; }));
    mesh::add_to_group(
      mesh,
      "pos_y",
      mesh::NodeSelection(
        [val = max_coord[1]](const mesh::NodeConstAccessor& n) -> bool { return n.y() == val; }));
    mesh::add_to_group(
      mesh,
      "neg_z",
      mesh::NodeSelection(
        [val = min_coord[2]](const mesh::NodeConstAccessor& n) -> bool { return n.z() == val; }));
    mesh::add_to_group(
      mesh,
      "pos_z",
      mesh::NodeSelection(
        [val = max_coord[2]](const mesh::NodeConstAccessor& n) -> bool { return n.z() == val; }));
}

static void add_unit_tension_x_bcs(BVPData& bvp_data)
{
    add_extent_node_groups(bvp_data.mesh());

    bvp_data.add_bc(DirichletBC("neg_x", { 0, 0 }));
    bvp_data.add_bc(DirichletBC("neg_y", { 1, 0 }));
    bvp_data.add_bc(DirichletBC("neg_z", { 2, 0 }));
    bvp_data.add_bc(DirichletBC("pos_x", { 0, 1 }));
}

static void add_unit_tension_z_bcs(BVPData& bvp_data)
{
    add_extent_node_groups(bvp_data.mesh());

    bvp_data.add_bc(DirichletBC("neg_x", { 0, 0 }));
    bvp_data.add_bc(DirichletBC("neg_y", { 1, 0 }));
    bvp_data.add_bc(DirichletBC("neg_z", { 2, 0 }));
    bvp_data.add_bc(DirichletBC("pos_z", { 2, 1 }));
}

static void add_unit_tension_x_z_bcs(BVPData& bvp_data)
{
    add_extent_node_groups(bvp_data.mesh());

    bvp_data.add_bc(DirichletBC("neg_x", { 0, 0 }));
    bvp_data.add_bc(DirichletBC("neg_y", { 1, 0 }));
    bvp_data.add_bc(DirichletBC("neg_z", { 2, 0 }));
    bvp_data.add_bc(DirichletBC("pos_x", { 0, 1 }));
    bvp_data.add_bc(DirichletBC("pos_z", { 2, 1 }));
}
}