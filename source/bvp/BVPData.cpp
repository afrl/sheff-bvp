#include "bvp/bvp/BVPData.h"

#include "bvp/bvp/Serialization.h"
#include "utility/MPIUtils.h"
#include "utility/serialization/Serialization.h"

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

namespace bvp
{
using json = nlohmann::json;

struct OffsetInformation
{
    int n_ranks;
    std::vector<std::streamoff> offsets;

    static void to_json(json& j, const OffsetInformation& info)
    {
        j = json { { "n_ranks", info.n_ranks }, { "offsets", info.offsets } };
    }

    static void from_json(const json& j, OffsetInformation& info)
    {
        j.at("n_ranks").get_to(info.n_ranks);
        j.at("offsets").get_to(info.offsets);
    }
};

BVPData load_bvp_data(const std::string filestem)
{
    return load_bvp_data(filestem, utility::MPIUtils::get_local_rank());
}

BVPData load_bvp_data(const std::string filestem, int rank_to_load)
{
    // Strip an extension if one was given
    fs::path filestem_fs    = fs::path(filestem).stem();
    std::string binary_path = filestem_fs.generic_string() + ".bvp.bin";
    std::string json_path   = filestem_fs.generic_string() + ".bvp.json";

    _UTIL_ASSERT(fs::exists(binary_path),
                 "The binary part of the data file \"" + binary_path + "\" did not exist.");
    _UTIL_ASSERT(fs::exists(json_path),
                 "The JSON part of the data file \"" + json_path + "\" did not exist.");

    json j;
    std::ifstream json_in(json_path);
    json_in >> j;
    OffsetInformation offset_info;
    OffsetInformation::from_json(j["offset_info"], offset_info);

    _UTIL_ASSERT(rank_to_load >= 0, "The rank to load must be greater than 0.");
    _UTIL_ASSERT(rank_to_load < offset_info.n_ranks,
                 "The rank to load was larger than the number of ranks that wrote the data file.");

    std::ifstream bin_in(binary_path, std::ios::in | std::ios::binary);
    utility::InputDeserializer s(bin_in);

    BVPData bvp_data;
    s(bvp_data);

    bin_in.close();

    return std::move(bvp_data);

}

void save_bvp_data(const BVPData& bvp_data, const std::string filestem)
{
    // Strip an extension if one was given
    fs::path filestem_fs    = fs::path(filestem).stem();
    std::string binary_path = filestem_fs.generic_string() + ".bvp.bin";
    std::string json_path   = filestem_fs.generic_string() + ".bvp.json";

    // Count how many bytes will be needed for the serialized data by rank
    utility::SerialByteCounter byte_counter;
    byte_counter(bvp_data);
    std::size_t n_bytes = byte_counter.n_bytes();
    auto size_by_rank   = utility::MPIUtils::gather_vectors<std::size_t>({ n_bytes });

    OffsetInformation offset_info;
    if (utility::MPIUtils::get_local_rank() == 0)
    {
        offset_info.n_ranks        = utility::MPIUtils::get_num_local_ranks();
        std::streamoff size_needed = 0;
        for (auto& s : size_by_rank)
        {
            offset_info.offsets.push_back(size_needed);
            size_needed += s;
        }

        // For now, zeroes will be written to a file to preallocate it.
        // the OS should have a way to preallocate a non-fragmented file
        // without writing the data, but it is not clear how to do it in a
        // portable way.  So for now...the dumb method will suffice.
        std::ofstream bin_stream(binary_path, std::ios::binary);
        const int buffer_size = 1024 * 1024; // 1 MB buffer
        std::vector<char> zero_buffer(buffer_size, 0);
        std::streamoff num_bytes_left = size_needed;
        while (num_bytes_left != 0)
        {
            if (num_bytes_left < buffer_size)
            {
                bin_stream.write(&zero_buffer[0], num_bytes_left);
                num_bytes_left = 0;
            }
            else
            {
                bin_stream.write(&zero_buffer[0], buffer_size);
                num_bytes_left -= buffer_size;
            }
        }
        bin_stream.close();

        // Write the JSON file for the meta information
        json j;
        OffsetInformation::to_json(j["offset_info"], offset_info);
        std::ofstream json_stream(json_path);
        json_stream << j;
        json_stream.close();
    }
    // Be sure the root rank finishes writing the file before opening it
    utility::MPIUtils::barrier();

    auto offset = utility::MPIUtils::scatter(offset_info.offsets);

    // Everyone needs to open the file in in/out mode to be able to seek to different parts and
    // write at will. Note: append mode will not let you go back before the end.
    std::ofstream out_stream;
    out_stream.open(binary_path, std::ios::binary | std::fstream::in | std::fstream::out);
    out_stream.seekp(offset, std::ios::beg);

    // Let each rank serialize their bvp_data
    utility::OutputSerializer s(out_stream);
    s(utility::named_ref("bvp_data", bvp_data));

    out_stream.close();

    utility::MPIUtils::barrier();
}

void save_bvp_data_json(const BVPData& bvp_data, const std::string filestem)
{
    auto filestem_fs = fs::path(filestem).stem();
    std::string filename;
    if (!utility::MPIUtils::is_distributed())
    {
        filename = filestem_fs.generic_string() + ".txt";
    }
    else
    {
        filename = filestem_fs.generic_string() + "_" +
                   std::to_string(utility::MPIUtils::get_local_rank()) + ".txt";
    }
    std::ofstream file(filename);
    utility::JSONSerializer s(file);
    s(utility::named_ref("bvp_data", bvp_data));
    s.write_out();
    file.close();
}

}