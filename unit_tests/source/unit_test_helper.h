#pragma once

#include <string>

#ifdef HAVE_FILESYSTEM
#include <experimental/filesystem>
#else
#include "fea_io/common/Filesystem.h"
#endif

namespace fea_analysis
{
    namespace fs = std::experimental::filesystem;

    /// Given the file name for a test data file, this returns the full filepath
	static const std::string get_test_data_filepath(
		const std::string& data_filename)
	{
        // Assumes test data is one directory up in another directory called "data"
		return fs::path(__FILE__).parent_path().parent_path().append("data").append(data_filename).string();
	}
}