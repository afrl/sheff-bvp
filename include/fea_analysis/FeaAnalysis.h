#pragma once

namespace fea_analysis
{
    class FeaAnalysis
    {
    public:
        FeaAnalysis() {}
        virtual ~FeaAnalysis() {}
        
        /// Loops over the time steps, does iteration, writes output, etc.
        /// It is the primary analysis method.
        void solve_time_steps();
        /// Restart from a previously output restart file.
        void restart_from_file(const std::string filename);
        
    protected:
        flex_mesh::Mesh<flex_mesh::Node, finite_element::FiniteElement> mesh_;
    
        /// Initialize the analysis 
        void initialize_analysis() = 0;
        void apply_boundary_conditions() = 0;
        void solve_time_step() = 0;
        void write_output() = 0;
        void write_restart_file();
    };
}