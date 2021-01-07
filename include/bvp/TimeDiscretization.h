#pragma once

#include "utility/Filesystem.h"

namespace bvp
{
enum class TimeDiscretizationType
{
    STATIC,
    LINEARLY_SPACED
};

/// Note: this struct needs to be expanded for psuedotime/time problems
struct TimeDiscretizationData
{
    int last_converged_time_step = -1;
    double last_convered_time    = 0.0;
    TimeDiscretizationType type  = TimeDiscretizationType::STATIC;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Serailization methods

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_ref("last_converged_time_step", last_converged_time_step),
          utility::named_ref("last_convered_time", last_convered_time),
          utility::named_value("type", type));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(last_converged_time_step, last_convered_time, type);
    }
};
}