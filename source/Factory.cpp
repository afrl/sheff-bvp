#include "bvp/Factory.h"

#include "bvp/bvp/FractureBVP.h"
#include "bvp/bvp/LinearElasticBVP.h"
#include "utility/StringUtils.h"

namespace bvp
{
BoundaryValueProblem& bvp_factory(const std::string& analysis_type,
                                  const std::string& kinematics_type)
{
    if (utility::istr_is_equal(analysis_type, "static_elasticity") &&
        utility::istr_is_equal(kinematics_type, "infinitesimal_strain"))
    {
        return LinearElasticBVP::get();
    }
    else if (utility::istr_is_equal(analysis_type, "static_elasticity_xfem") &&
             utility::istr_is_equal(kinematics_type, "infinitesimal_strain"))
    {
        return FractureBVP::get();
    }
    throw _UTIL_MSG_EXCEPT("Could not create an analysis of type " + analysis_type +
                           " with kinematics " + kinematics_type);
}
}