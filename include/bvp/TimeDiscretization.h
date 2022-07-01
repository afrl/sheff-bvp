#pragma once

#include "simple_math/EigenDefines.h"
#include "utility/Filesystem.h"

namespace bvp
{
enum class ScaleFunctionType
{
    STATIC,
    LINEARLY_SPACED
};

/// Notes:
/// * For any time less than times.front() or any time greater than times.back(), the scale factor
///   will be assumed to be zero.
/// * The length of scale_factors and times must match, since the are effectively (x,y) pairs that
///   define the scale function.
/// * For a given time between the values in "times", the scale factor will be interpolated assuming
///   a linear variation between the two respective values in "scale_factors".
struct ScaleFunction
{
    /// An identifier to use in the input file
    std::string name;
    // The x-axis of the scale function
    std::vector<double> times;
    // The y-axis of the scale function
    std::vector<double> scale_factors;

    double get_scale_factor(double time)
    {
        _UTIL_ASSERT(times.size() > 1,
                     "There must be at least one time in the scale function: " + name);
        if (time < times.front() || time > times.back())
        {
            return 0.0;
        }
        else
        {
            for (auto i : utility::range(times.size() - 1))
            {
                if (time >= times[i])
                {
                    const double interval    = times[i + 1] - times[i];
                    const double weight_i    = (times[i + 1] - time) / interval;
                    const double weight_i_p1 = (time - times[i]) / interval;
                    return weight_i * scale_factors[i] + weight_i_p1 * scale_factors[i + 1];
                }
            }
            // If we made it hear, then this is the only possibility
            _UTIL_ASSERT(times.back() == time,
                         "Something is wrong with the logic of this function.");
            return scale_factors.back();
        }
    }
};

/// Note: this struct needs to be expanded for psuedotime/time problems
struct TimeDiscretizationData
{
    /// The last time step that was converged
    int last_converged_time_step = -1;
    double last_converged_time   = 0.0;
    double current_delta_time    = 1.0;
    double current_time          = 0.0;

    /// The list of times given by the user that should be solved.
    std::vector<double> target_times = { 1.0 };

    /// The list of scale functions that can be used for boundary conditions
    std::vector<ScaleFunction> scale_functions;

    // Substepping variables

    /// Number of substeps that have converged without a cut back
    int n_consecutive_converged_substeps = 0;
    /// Number of substeps to converge without a cut back before larger increments are taken
    int n_consecutive_substeps_before_increase = 0;
    /// Factor to cutback time increment if substep fails to converge (delta_time /= cutback_factor)
    int cutback_factor = 2;
    /// Factor to increase time increment if n_consecutive_converged_substeps ==
    /// n_consecutive_substeps_before_increase (delta_time *= cutback_factor)
    int accelerate_factor = 2;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Serailization methods

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_ref("last_converged_time_step", last_converged_time_step),
          utility::named_ref("last_converged_time", last_converged_time),
          utility::named_ref("current_delta_time", current_delta_time),
          utility::named_ref("current_time", current_time),
          /* TODO write target_times */
          /* TODO write scale_functions */
          utility::named_ref("n_consecutive_converged_substeps", n_consecutive_converged_substeps),
          utility::named_ref("n_consecutive_substeps_before_increase",
                             n_consecutive_substeps_before_increase),
          utility::named_ref("cutback_factor", cutback_factor),
          utility::named_ref("accelerate_factor", accelerate_factor));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(last_converged_time_step,
          last_converged_time,
          current_delta_time,
          current_time,
          /* TODO write target_times */
          /* TODO write scale_functions */
          n_consecutive_converged_substeps,
          n_consecutive_substeps_before_increase,
          cutback_factor,
          accelerate_factor);
    }
};
}