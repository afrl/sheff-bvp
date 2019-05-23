#include "fea_analysis/LinearAnalysis.h"

#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/Vector.h"
#include "sparse_linear_algebra/SparseSolver.h"

#include "fea_crack/RXFEMCrack.h"

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

        /// Creates all of the objects needed before solve_time_steps can be called
        virtual void setup_analysis() override;

        void insert_precracks(const std::vector<simple_math::Point>& initiation_points);

    protected:
        fea_crack::RXFEMCrack::ElementAdjacencyList adjacency_list_;
        std::vector<std::unique_ptr<fea_crack::RXFEMCrack>> cracks_;

        virtual void solve_time_step() override;

        void solve_equillibrium_step();
        bool try_to_initiate_cracks();
        bool try_to_grow_cracks();
    };
}