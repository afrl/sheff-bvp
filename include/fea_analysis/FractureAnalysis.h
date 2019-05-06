#include "fea_analysis/LinearAnalysis.h"

#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/Vector.h"
#include "sparse_linear_algebra/SparseSolver.h"

#include <memory>
#include <vector>

namespace fea_analysis
{
    /// A linear elastic analysis
    class FractureAnalysis : public LinearElasticAnalysis
    {
    public:
        FractureAnalysis() {}
        virtual ~FractureAnalysis() {}

    protected:
        virtual void solve_time_step() override;
    };
}