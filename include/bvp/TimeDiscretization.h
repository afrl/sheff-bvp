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
    // Stores the last scale factor for the last converged substep in order to calculate the
    // incremental scale factor (current - last value)
    double scale_factor_for_last_converged_substep = 0.0;

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


    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_ref("name", name),
          utility::named_value("n_times", times.size()),
          utility::named_value("times", utility::slice_c(times)),
          utility::named_value("n_scale_factors", scale_factors.size()),
          utility::named_value("scale_factors", utility::slice_c(scale_factors)),
          utility::named_ref("scale_factor_for_last_converged_substep",
                             scale_factor_for_last_converged_substep));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        std::size_t n_times, n_scale_factors;
        s(name, n_times);
        times.resize(n_times);
        s(utility::slice(times), n_scale_factors);
        scale_factors.resize(n_scale_factors);
        s(utility::slice(scale_factors), scale_factor_for_last_converged_substep);
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
          utility::named_value("n_target_times", target_times.size()),
          utility::named_value("target_times", utility::slice_c(target_times)),
          utility::named_value("n_scale_functions", scale_functions.size()),
          utility::named_value("scale_functions", utility::slice_c(scale_functions)),
          utility::named_ref("n_consecutive_converged_substeps", n_consecutive_converged_substeps),
          utility::named_ref("n_consecutive_substeps_before_increase",
                             n_consecutive_substeps_before_increase),
          utility::named_ref("cutback_factor", cutback_factor),
          utility::named_ref("accelerate_factor", accelerate_factor));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(last_converged_time_step, last_converged_time, current_delta_time, current_time);
        
        std::size_t n_target_times;
        s(n_target_times);
        target_times.resize(n_target_times);
        s(utility::slice(target_times));

        std::size_t n_scale_functions;
        s(n_scale_functions);
        scale_functions.resize(n_scale_functions);
        s(utility::slice(scale_functions));
        
        s(n_consecutive_converged_substeps,
          n_consecutive_substeps_before_increase,
          cutback_factor,
          accelerate_factor);
    }
};
}