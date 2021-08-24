#include "bvp/bvp/BoundaryValueProblem.h"

#include "assembler/SerialAssembler.h"
#include "bvp/bvp/Serialization.h"
#include "constitutive/field/FieldProperties.h"
#include "fea_io/meta/MetaReader.h"
#include "fea_io/meta/MetaWriter.h"
#include "finite_element/ModelSelectionHelper.h"
#include "mesh/algorithms/IO.h"
#include "utility/MPIUtils.h"
#include "utility/Profile.h"
#include "utility/serialization/Serialization.h"

namespace bvp
{
static void copy_file(fs::path old_path, fs::path new_path, const std::vector<std::string>& exts)
{
    if (utility::MPIUtils::get_local_rank() == 0)
    {
        for (auto& ext : exts)
        {
            auto old_path_ext = old_path;
            old_path_ext += ext;
            auto new_path_ext = new_path;
            new_path_ext += ext;
            fs::copy_file(old_path_ext, new_path_ext, fs::copy_options::overwrite_existing);
        }
    }
}

void BoundaryValueProblem::set_input(const std::string input_script_path,
                                     const std::string input_meta_file_path,
                                     const std::string output_dir_path,
                                     BVPData& bvp_data) const
{
    fea_io::MetaReader reader;
    reader.open(input_meta_file_path);
    auto& input_analysis_meta = reader.get_analysis_meta();

    fs::path input_meta_dir = fs::path(input_meta_file_path).parent_path();
    fs::path output_dir     = output_dir_path;
    if (!fs::exists(output_dir))
    {
        fs::create_directory(output_dir);
    }
    else if (!fs::is_directory(output_dir))
    {
        throw _UTIL_MSG_EXCEPT_C(
          "The path given for the output directory exists but is not a directory.");
    }
    bvp_data.output_dir() = output_dir;

    // Copy over any mesh, group, and data files into the output directory in order to keep
    // the directory self-contained.  While the files are copied over, populate the InputInfo
    // struct for the new output meta file with information from the input meta file.
    fea_io::InputInfo output_input_info;      // InputInfo for the output meta file
    output_input_info.input_script_path = ""; // TODO
    for (auto& input_lod : input_analysis_meta.input_info.levels_of_detail)
    {
        fea_io::LevelOfDetail output_lod;
        output_lod = input_lod;
        output_lod.mesh_objects.clear();

        for (auto& mesh_object : input_lod.mesh_objects)
        {
            fea_io::MeshObject output_mesh_object;

            fs::path mesh_path = input_meta_dir / mesh_object.mesh_path;
            _UTIL_ASSERT_C(mesh_path.parent_path() == input_meta_dir,
                           "All paths in the input meta file must not point to a directory outside "
                           "where the meta file is located.");
            copy_file(mesh_path, output_dir / mesh_path.filename(), { ".m.json", ".m.bin" });
            output_mesh_object.mesh_path = mesh_path.filename().string();

            for (auto& data_path_str : mesh_object.data_paths)
            {
                fs::path data_path = input_meta_dir / data_path_str;
                _UTIL_ASSERT_C(data_path.parent_path() == input_meta_dir,
                               "All paths in the input meta file must not point to a directory "
                               "outside where the meta file is located.");
                copy_file(data_path, output_dir / data_path.filename(), { ".d.json", ".d.bin" });
                output_mesh_object.data_paths.push_back(data_path.filename().string());
            }

            for (auto& group_path_str : mesh_object.group_paths)
            {
                fs::path group_path = input_meta_dir / group_path_str;
                _UTIL_ASSERT_C(group_path.parent_path() == input_meta_dir,
                               "All paths in the input meta file must not point to a directory "
                               "outside where the meta file is located.");
                copy_file(group_path, output_dir / group_path.filename(), { ".g.json", ".g.bin" });
                output_mesh_object.group_paths.push_back(group_path.filename().string());
            }

            output_mesh_object.type = mesh_object.type;

            output_lod.mesh_objects.push_back(output_mesh_object);
        }

        output_input_info.levels_of_detail.push_back(output_lod);
    }

    // Read the mesh into the mesh data structure for the analysis
    _UTIL_ASSERT_C(input_analysis_meta.input_info.levels_of_detail.size() == 1,
                   "Currently, only one level of detail is supported for the input data.  Relaxing "
                   "this restriction will require significant development.");
    _UTIL_ASSERT_C(input_analysis_meta.input_info.levels_of_detail[0].mesh_objects.size() == 1,
                   "Currently, only one mesh object is supported for the input data.  Relaxing "
                   "this restriction will require development.");

    std::string mesh_filepath =
      (input_meta_dir /
       input_analysis_meta.input_info.levels_of_detail[0].mesh_objects[0].mesh_path)
        .string();

    auto mesh = mesh::read_fea_io_from_meta(mesh_filepath, utility::MPIUtils::get_local_rank());
    bvp_data.fe_mesh().initialize(std::move(mesh));

    // For debugging TODO delete
    // flex_mesh::write_vtk(*mesh_, "partition_" +
    // std::to_string(utility::MPIUtils::get_local_rank()) + ".vtk", true, false);

    fea_io::MetaWriter output_meta;
    output_meta.open((bvp_data.output_dir() / "output").string());
    auto& output_analysis_meta      = output_meta.get_analysis_meta();
    output_analysis_meta            = input_analysis_meta;
    output_analysis_meta.input_info = output_input_info;

    output_meta.update_meta_file();
}

void BoundaryValueProblem::set_input(mesh::Mesh&& mesh,
                                     const std::string analysis_name,
                                     const std::string output_dir_path,
                                     BVPData& bvp_data) const
{
    bvp_data.fe_mesh().initialize(std::move(mesh));

    fs::path output_dir = output_dir_path;
    if (!fs::exists(output_dir))
    {
        fs::create_directory(output_dir);
    }
    else if (!fs::is_directory(output_dir))
    {
        throw _UTIL_MSG_EXCEPT_C(
          "The path given for the output directory exists but is not a directory.");
    }
    bvp_data.output_dir() = output_dir;

    mesh::write_fea_io(bvp_data.mesh(), output_dir.generic_string(), "input_mesh");

    fea_io::MetaWriter output_meta;
    output_meta.open((output_dir / "output").string());
    auto& output_analysis_meta = output_meta.get_analysis_meta();
    output_analysis_meta.name  = analysis_name;
    output_analysis_meta.levels_of_detail_names.push_back("original");

    fea_io::MeshObject mesh_object;
    mesh_object.mesh_path = "input_mesh";
    mesh_object.type      = fea_io::MeshObjectType::Volume;

    output_analysis_meta.input_info.levels_of_detail.push_back(
      { 0, fea_io::LODParentRelation::None, { mesh_object } });

    output_meta.update_meta_file();
}

void BoundaryValueProblem::set_field_models(finite_element::BasalFieldSelection selection,
                                            const std::string& formulation_name,
                                            const std::string& properties_name,
                                            BVPData& bvp_data) const
{
    selection.run_selection(bvp_data.fe_mesh());
    std::vector<finite_element::SelectedBasalField> selected_fields(
      selection.get_selection().begin(),
      selection.get_selection().end());

    set_field_models(selected_fields, formulation_name, properties_name, bvp_data);
}

void BoundaryValueProblem::set_field_models(
  const std::vector<finite_element::SelectedBasalField>& selected_fields,
  const std::string& formulation_name,
  const std::string& properties_name,
  BVPData& bvp_data) const
{
    _F_BEGIN();

    if (bvp_data.mesh().elements().size() == 0)
    {
        throw _UTIL_MSG_EXCEPT_C("The mesh must be set before models can be assigned.");
    }

    auto props_ref = constitutive::FieldPropertiesCollection::singleton().get_ref(properties_name);

    auto& selection_helper = finite_element::ModelSelectionHelper::get();
    for (auto selected_field : selected_fields)
    {
        selection_helper.assign_field_model(
          bvp_data.fe_mesh().get_finite_element(selected_field.element_index),
          selected_field.basal_field_index,
          formulation_name,
          props_ref);
    }

    _F_END();
}

void BoundaryValueProblem::solve_time_steps(BVPData& bvp_data) const
{
    _F_BEGIN2(bvp_data.bvp_name());

    init_all_data(bvp_data);

    // Time loop
    SolveTimeStepResult solve_result;
    while (!solve_result.end_of_time_steps)
    {
        solve_result = solve_next_time_step(bvp_data);

        if (!solve_result.converged)
        {
            _UTIL_WARN_C("Time step failed to converge, exiting solve_time_steps");
            break;
        }
    }

    if (sparse_linear_algebra::LinAlgFactory::get().verbose_level() >=
        sparse_linear_algebra::VerboseLevel::BASIC)
    {
        sparse_linear_algebra::LinAlgFactory::get().write_log(
          (bvp_data.output_dir() / "lin_alg_log.txt").generic_string());
    }

    _F_END();
}

void BoundaryValueProblem::init_core_data(BVPData& bvp_data) const
{
    _F_BEGIN();

    // Setup the dof map and any global matrices/vectors.  The dof map will be
    // used to setup the containers.
    if (bvp_data.dof_manager().get_dof_maps().n_dofs_on_rank() == 0)
    {
        if (!bvp_data.adjacency_list().initialized())
        {
            bvp_data.node_to_element_map().update(bvp_data.mesh());
            bvp_data.adjacency_list().update(bvp_data.mesh(), bvp_data.node_to_element_map());
        }
        if (!bvp_data.element_scatter_map().initialized())
        {
            bvp_data.element_scatter_map().update(bvp_data.fe_mesh(), bvp_data.adjacency_list());
        }

        bvp_data.dof_manager().enumerate(bvp_data.fe_mesh(),
                                         bvp_data.adjacency_list(),
                                         bvp_data.element_scatter_map());
    }

    if (bvp_data.core_data.solution.empty())
    {
        bvp_data.core_data.solution.resize(bvp_data.dof_manager().get_dof_maps().n_dofs_on_rank(),
                                           0);
    }

    // Note: allocate_isvs will only allocate if the field does not already have ISV's allocated.
    for (auto fe_accessor : bvp_data.fe_mesh().finite_elements())
    {
        finite_element::FiniteElement::allocate_isvs(fe_accessor.layout);
    }

    _F_END();
}

void BoundaryValueProblem::init_all_data(BVPData& bvp_data) const
{
    _F_BEGIN();

    init_core_data(bvp_data);

    if (!bvp_data.adjacency_list().initialized())
    {
        bvp_data.node_to_element_map().update(bvp_data.mesh());
        bvp_data.adjacency_list().update(bvp_data.mesh());
    }
    if (!bvp_data.element_scatter_map().initialized())
    {
        bvp_data.element_scatter_map().update(bvp_data.fe_mesh(), bvp_data.adjacency_list());
    }
    if (!bvp_data.secondary_data.assembler_ptr)
    {
        bvp_data.secondary_data.assembler_ptr =
          std::make_unique<assembler::SerialAssembler>(bvp_data.fe_mesh());
    }

    if (!bvp_data.gl_lin_alg)
    {
        init_lin_algebra(bvp_data);
    }

    _F_END();
}
}