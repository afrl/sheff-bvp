#pragma once

#include "fea_analysis/FeaAnalysis.h"

namespace fea_analysis
{
	/// Creates the given type of analysis.  Supports:
	/// - "linear_elastic"
	/// - "fracture"
	std::unique_ptr<FeaAnalysis> analysis_factory(const std::string& type);
}