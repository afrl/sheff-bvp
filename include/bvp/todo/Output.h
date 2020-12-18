#pragma once

#include "output/Output.h"

namespace bvp
{
struct OutputData
{
    std::vector<std::unique_ptr<output::OutputFlag<finite_element::VolumeQuantity>>>
      volume_output_flags;
    std::vector<std::unique_ptr<output::OutputFlag<finite_element::SurfaceQuantity>>>
      surface_output_flags;
    fea_io::MetaWriter output_meta;
};
}