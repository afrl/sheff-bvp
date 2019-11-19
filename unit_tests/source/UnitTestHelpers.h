#pragma once

#include "utility/Filesystem.h"
#include "fea_analysis/BoundaryCondition.h"

namespace fea_analysis
{
    /// Given the file name for a test data file, this returns the full filepath
    static const std::string get_test_data_filepath(
        const std::string& data_filename)
    {
        // Assumes test data is one directory up in another directory called "data"
        auto data_path = fs::path(__FILE__).parent_path().parent_path().append("data").append(data_filename);
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
                throw _UTIL_MSG_EXCEPT("The directory for output (" + output_dir.string() + ") already exists and is not a directory.");
            }
            else
            {
                fs::remove_all(output_dir);
                fs::create_directory(output_dir);
            }
        }
    }

    static std::vector<std::unique_ptr<DirichletBC>> get_unit_tension_x_bcs()
    {
        std::vector<std::unique_ptr<DirichletBC>> bcs;

        // Fix u along the x = -1 side
        NodeSelection x_neg_1_selection;
        x_neg_1_selection.add_criteria(
            [](const fea_analysis::Node& n) -> bool {
            return n.x() == -1.0;
        });
        bcs.push_back(std::unique_ptr<DirichletBC>(
            new DirichletBC(x_neg_1_selection, { 0 , 0 })));
        // Fix v along the y = -1 side
        NodeSelection y_neg_1_selection;
        y_neg_1_selection.add_criteria(
            [](const fea_analysis::Node& n) -> bool {
            return n.y() == -1.0;
        });
        bcs.push_back(std::unique_ptr<DirichletBC>(
            new DirichletBC(y_neg_1_selection, { 1 , 0 })));
        // Fix w along the z = -1 side
        NodeSelection z_neg_1_selection;
        z_neg_1_selection.add_criteria(
            [](const fea_analysis::Node& n) -> bool {
            return n.z() == -1.0;
        });
        bcs.push_back(std::unique_ptr<DirichletBC>(
            new DirichletBC(z_neg_1_selection, { 2 , 0 })));
        // Specify u along the x = 1 side
        NodeSelection x_pos_1_selection;
        x_pos_1_selection.add_criteria(
            [](const fea_analysis::Node& n) -> bool {
            return n.x() == 1.0;
        });
        bcs.push_back(std::unique_ptr<DirichletBC>(
            new DirichletBC(x_pos_1_selection, { 0 , 1 })));

        return bcs;
    }
}