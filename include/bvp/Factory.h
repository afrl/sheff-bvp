#pragma once

#include "bvp/bvp/BoundaryValueProblem.h"

namespace bvp
{
/// Creates the given type of analysis.
/// \param analysis_type - the type of analysis.  Currently supports:
/// - "static_elasticity"
/// - "static_elasticity_xfem"
/// \param kinematics_type - the type of kinematics for the analysis (if applicable)
BoundaryValueProblem& bvp_factory(const std::string& analysis_type,
                                  const std::string& kinematics_type);
}