#include "fea_analysis/FeaAnalysis.h"

namespace fea_analysis
{
    void FeaAnalysis::set_mesh(std::unique_ptr<MeshType> mesh)
    {
        mesh_ = std::move(mesh);
    }

    void FeaAnalysis::add_constitutive_relation(const std::string name, std::unique_ptr<constitutive::BaseConstitutiveRelation> relation)
    {
        constitutive_relations_[name] = std::move(relation);
    }

    void FeaAnalysis::solve_time_steps()
    {
        // If no time step manager has been set through the input script, it is assumed that
        // this is a static analysis.
        if (time_step_manager_ == 0)
            time_step_manager_ = std::make_unique<StaticAnalysisTimeStepManager>();

        // Setup the dof map and any global matrices/vectors.  The dof map will be
        // used to setup the containers.
        dof_map_ = std::make_unique<dof_manager::DofMap>(*mesh_);
        initialize_linear_algebra();

        // Time loop
        time_step_manager_->start_time_step_loop();
        do
        {
            /* TODO add restart capability */
            solve_time_step();
            write_output();
        } while (time_step_manager_->next_time_step());
    }
}