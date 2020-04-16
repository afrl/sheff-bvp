#include "fea_analysis/Factory.h"

#include "fea_analysis/LinearAnalysis.h"
#include "fea_analysis/FractureAnalysis.h"

#include "utility/StringUtils.h"

namespace fea_analysis
{
	std::unique_ptr<FeaAnalysis> analysis_factory(
		const std::string& analysis_type,
		const std::string& kinematics_type)
	{
		if (utility::istr_is_equal(analysis_type, "static_elasticity") &&
			utility::istr_is_equal(kinematics_type, "infinitesimal_strain"))
		{
			return std::unique_ptr<FeaAnalysis>(new LinearElasticAnalysis());
		}
		else if (utility::istr_is_equal(analysis_type, "fracture") &&
				 utility::istr_is_equal(kinematics_type, "infinitesimal_strain"))
		{
			return std::unique_ptr<FeaAnalysis>(new FractureAnalysis());
		}
		throw _UTIL_MSG_EXCEPT("Could not create an analysis of type " + analysis_type + " with kinematics " + kinematics_type);
	}
}