#include "fea_analysis/FeaAnalysis.h"

#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/Vector.h"
#include "sparse_linear_algebra/SparseSolver.h"

namespace fea_analysis
{
    /// A linear elastic analysis
    class LinearElasticAnalysis : public FeaAnalysis
    {
    public:
        LinearElasticAnalysis() {}
        virtual ~LinearElasticAnalysis() {}

    protected:
        std::unique_ptr<sparse_linear_algebra::SparseMatrix> gl_stiffness_;
        std::unique_ptr<sparse_linear_algebra::Vector> gl_displacement_;
        std::unique_ptr<sparse_linear_algebra::Vector> gl_eigenstrain_force_;
        std::unique_ptr<sparse_linear_algebra::Vector> gl_rhs_;
        std::unique_ptr<sparse_linear_algebra::SparseSolver> sparse_solver_;

        /// Initialize the global matrices and vectors for the analysis
		virtual void set_element_formulations() override;
        virtual void initialize_linear_algebra() override;
        virtual void apply_boundary_conditions() override {}
        virtual void solve_time_step() override;
        virtual void write_output() override {}
        virtual void write_restart_file() override {}
    };
}