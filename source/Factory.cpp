#include "fea_analysis/Factory.h"

#include "fea_analysis/LinearAnalysis.h"
#include "fea_analysis/FractureAnalysis.h"

#include "utility/StringUtils.h"

namespace fea_analysis
{
	std::unique_ptr<FeaAnalysis> analysis_factory(const std::string& type)
	{
		if (utility::istr_is_equal(type, "linear_elastic"))
		{
			return std::unique_ptr<FeaAnalysis>(new LinearElasticAnalysis());
		}
		else if (utility::istr_is_equal(type, "fracture"))
		{
			return std::unique_ptr<FeaAnalysis>(new FractureAnalysis());
		}
		throw _UTIL_MSG_EXCEPT("Could not create an analysis of type: " + type);
	}
}