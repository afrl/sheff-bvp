#pragma once

#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "finite_element/FiniteElement.h"
#include "constitutive/BaseConstitutiveRelation.h"
#include "assembler/Assembler.h"
#include "dof_manager/DofMap.h"

#include "TimeStepManager.h"

namespace fea_analysis
{
    class FeaAnalysis
    {
    public:
        using MeshType = flex_mesh::Mesh<flex_mesh::Node, finite_element::FiniteElement>;

        FeaAnalysis() {}
        virtual ~FeaAnalysis() {}

        /// Sets the mesh for the analysis
        void set_mesh(std::unique_ptr<MeshType> mesh);
        /// Add a constitutive relation to the analysis
        void add_constitutive_relation(std::string, std::unique_ptr<constitutive::BaseConstitutiveRelation>);

        /// Loops over the time steps, does iteration, writes output, etc.
        /// It is the primary analysis method.
        void solve_time_steps();
        /// Restart from a previously output restart file.
        void restart_from_file(const std::string filename);

    protected:
        std::unique_ptr<TimeStepManager> time_step_manager_;
        std::unique_ptr<MeshType> mesh_;
        std::unique_ptr<assembler::Assembler> assembler_;
        std::unique_ptr<dof_manager::DofMap> dof_map_;
        //std::vector<std::unique_ptr<BoundaryCondition>> boundary_conditions_;
        std::map<std::string, std::unique_ptr<constitutive::BaseConstitutiveRelation>> constitutive_relations_;

        /// Initialize the global matrices and vectors for the analysis
        virtual void initialize_linear_algebra() = 0;
        virtual void apply_boundary_conditions() = 0;
        virtual void solve_time_step() = 0;
        virtual void write_output() = 0;
        virtual void write_restart_file() = 0;

        template<typename GlobalContainer>
        void initialize_container(GlobalContainer& gl_container);
    };
    template<typename GlobalContainer>
    inline void FeaAnalysis::initialize_container(GlobalContainer & gl_container)
    {
        gl_container.set_local_global_dof_map(
            dof_map_->get_element_to_global_map(),
            dof_map_->get_element_num_dofs(),
            mesh_->elements().size(),
            dof_map_->get_n_local_dofs(),
            dof_map_->get_n_global_dofs());
    }
}