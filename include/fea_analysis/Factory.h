#pragma once

#include "fea_analysis/FeaAnalysis.h"

namespace fea_analysis
{
	/// Creates the given type of analysis.
	/// \param analysis_type - the type of analysis.  Currently supports:
	/// - "static_elasticity"
	/// - "static_elasticity_xfem"
	/// \param kinematics_type - the type of kinematics for the analysis (if applicable)
	std::unique_ptr<FeaAnalysis> analysis_factory(
		const std::string& analysis_type,
		const std::string& kinematics_type);
}