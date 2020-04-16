#include "fea_analysis/FeaAnalysis.h"

#include "assembler/SerialAssembler.h"
#include "fea_io/meta/MetaReader.h"
#include "flex_mesh/algorithms/IOFeaIo.h"

#include "finite_element/Factory.h"
#include "constitutive/field/FieldProperties.h"

namespace fea_analysis
{
    static void copy_file(
        fs::path old_path,
        fs::path new_path,
        const std::vector<std::string>& exts)
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

    void FeaAnalysis::set_input(
        const std::string input_script_path,
        const std::string input_meta_file_path,
        const std::string output_dir_path)
    {
        fea_io::MetaReader reader;
        reader.open(input_meta_file_path);
        auto& input_analysis_meta = reader.get_analysis_meta();

        fs::path input_meta_dir = fs::path(input_meta_file_path).parent_path();
        fs::path output_dir = output_dir_path;
        if (!fs::exists(output_dir))
        {
            fs::create_directory(output_dir);
        }
        else if (!fs::is_directory(output_dir))
        {
            throw _UTIL_MSG_EXCEPT_C("The path given for the output directory exists but is not a directory.");
        }
        output_dir_ = output_dir;

        // Copy over any mesh, group, and data files into the output directory in order to keep
        // the directory self-contained.  While the files are copied over, populate the InputInfo
        // struct for the new output meta file with information from the input meta file.
        fea_io::InputInfo output_input_info; // InputInfo for the output meta file
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
                    "All paths in the input meta file must not point to a directory outside where the meta file is located.");
                copy_file(mesh_path, output_dir / mesh_path.filename(),
                    { ".m.json", ".m.bin" });
                output_mesh_object.mesh_path = mesh_path.filename().string();

                for (auto& data_path_str : mesh_object.data_paths)
                {
                    fs::path data_path = input_meta_dir / data_path_str;
                    _UTIL_ASSERT_C(data_path.parent_path() == input_meta_dir,
                        "All paths in the input meta file must not point to a directory outside where the meta file is located.");
                    copy_file(data_path, output_dir / data_path.filename(),
                        { ".d.json", ".d.bin" });
                    output_mesh_object.data_paths.push_back(data_path.filename().string());
                }

                for (auto& group_path_str : mesh_object.group_paths)
                {
                    fs::path group_path = input_meta_dir / group_path_str;
                    _UTIL_ASSERT_C(group_path.parent_path() == input_meta_dir,
                        "All paths in the input meta file must not point to a directory outside where the meta file is located.");
                    copy_file(group_path, output_dir / group_path.filename(),
                        { ".d.json", ".d.bin" });
                    output_mesh_object.data_paths.push_back(group_path.filename().string());
                }

                output_lod.mesh_objects.push_back(output_mesh_object);
            }

            output_input_info.levels_of_detail.push_back(output_lod);
        }
        
        // Read the mesh into the mesh data structure for the analysis
        _UTIL_ASSERT_C(input_analysis_meta.input_info.levels_of_detail.size() == 1,
            "Currently, only one level of detail is supported for the input data.  Relaxing this restriction will require significant development.");
        _UTIL_ASSERT_C(input_analysis_meta.input_info.levels_of_detail[0].mesh_objects.size() == 1,
            "Currently, only one mesh object is supported for the input data.  Relaxing this restriction will require development.");

        std::string mesh_filepath = (input_meta_dir /
            input_analysis_meta.input_info.levels_of_detail[0].mesh_objects[0].mesh_path).string();
        mesh_ = std::make_unique<Mesh>();
        flex_mesh::FEAIOInterface<Mesh>::read_partition(*mesh_, mesh_filepath, 0, 0);

        output_meta_.open("output");
        auto& output_analysis_meta = output_meta_.get_analysis_meta();
        output_analysis_meta = input_analysis_meta;
        output_analysis_meta.input_info = output_input_info;

        output_meta_.update_meta_file();
    }

    void FeaAnalysis::set_input(
        std::unique_ptr<Mesh> mesh,
        const std::string analysis_name,
        const std::string output_dir_path)
    {
        mesh_ = std::move(mesh);

        fs::path output_dir = output_dir_path;
        if (!fs::exists(output_dir))
        {
            fs::create_directory(output_dir);
        }
        else if (!fs::is_directory(output_dir))
        {
            throw _UTIL_MSG_EXCEPT_C("The path given for the output directory exists but is not a directory.");
        }
        output_dir_ = output_dir;

        std::string mesh_file_path = (output_dir / "input_mesh").string();
        if (utility::MPIUtils::is_distributed())
        {
            std::cout << "TODO: write mesh" << std::endl;
        }
        else
        {
            flex_mesh::FEAIOInterface<Mesh>::write(*mesh_, mesh_file_path);
        }

        output_meta_.open((output_dir_ / "output").string());
        auto& output_analysis_meta = output_meta_.get_analysis_meta();
        output_analysis_meta.name = analysis_name;
        output_analysis_meta.levels_of_detail_names.push_back("original");

        fea_io::MeshObject mesh_object;
        mesh_object.mesh_path = "input_mesh";
        mesh_object.type = fea_io::MeshObjectType::Volume;

        output_analysis_meta.input_info.levels_of_detail.push_back({
            0,
            fea_io::LODParentRelation::None,
            { mesh_object } });

        output_meta_.update_meta_file();
    }

    void FeaAnalysis::set_field_model(
        finite_element::BasalFieldSelection selection,
        const std::string& formulation_name,
        const constitutive::FieldProperties& properties)
    {
        if (mesh_ == nullptr)
        {
            throw _UTIL_MSG_EXCEPT_C("The mesh must be set before models can be assigned.");
        }

        field_properties_.push_back(properties.clone());

        selection.run_selection(*mesh_);

        auto& factory = finite_element::Factory::get();
        for (auto selected_field : selection.get_selection())
        {
            factory.assign_field_model(
                selected_field.element,
                selected_field.index,
                formulation_name,
                *field_properties_.back());
        }
    }

	void FeaAnalysis::add_boundary_condition(std::unique_ptr<BoundaryCondition> boundary_condition)
	{
		boundary_conditions_.push_back(std::move(boundary_condition));
	}

    void fea_analysis::FeaAnalysis::setup_analysis()
    {
        // If no time step manager has been set through the input script, it is assumed that
        // this is a static analysis.
        if (time_step_manager_ == 0)
            time_step_manager_ = std::make_unique<StaticAnalysisTimeStepManager>();

        assembler_ = std::make_unique<assembler::SerialAssembler>(*mesh_);

        // Setup the dof map and any global matrices/vectors.  The dof map will be
        // used to setup the containers.
        dof_map_ = std::make_unique<dof_manager::DofMap>(*mesh_);
        initialize_linear_algebra();

        for (auto& e_ptr : mesh_->elements())
        {
            e_ptr->allocate_isvs();
        }
    }

    void FeaAnalysis::solve_time_steps()
    {
        pre_solve_tasks();

        // Time loop
        time_step_manager_->start_time_step_loop();
        do
        {
            // Create a directory for the current time step
            auto curr_time_step = time_step_manager_->get_time_step();
            curr_time_step_dir_ = output_dir_ / ("ts_" + std::to_string(curr_time_step));
            fs::create_directory(curr_time_step_dir_);

            /* TODO add restart capability */
            solve_time_step();
            write_output();
        } while (time_step_manager_->next_time_step());
    }

    void FeaAnalysis::apply_boundary_conditions()
    {
        for (auto& bc : boundary_conditions_)
        {
            // IMPORTANT NOTE: this assumes the field has # of components as dimensionality
            // This works for displacement with solid elements, but reconsider what should
            // be done to handle scalar fields, shell elements, combined fields, etc.
            // TODO figure out a robust term strategy for this
            bc->initialize(*mesh_, *dof_map_,
                utility::GlobalVariables::get_dimensionality(),
                *gl_constraint_matrix_);
        }
    }
}