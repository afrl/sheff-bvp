#pragma once

#include "finite_element/FiniteElement.h"
#include "constitutive/BaseConstitutiveRelation.h"
#include "assembler/Assembler.h"
#include "dof_manager/DofMap.h"
#include "flex_mesh/maps/SubcomponentMaps.h"

#include "output/Output.h"
#include "fea_io/meta/MetaWriter.h"

#include "Defines.h"
#include "BoundaryCondition.h"
#include "TimeStepManager.h"

namespace fea_analysis
{
    class FeaAnalysis
    {
    public:
        FeaAnalysis() {}
        virtual ~FeaAnalysis() {}

        /// Read an input meta file (will set the mesh)
        void set_input_meta(const std::string relative_path_to_meta);
        /// Sets the mesh for the analysis
        void set_mesh(std::unique_ptr<Mesh> mesh);
        /// Add a constitutive relation to the analysis
        void add_constitutive_relation(const std::string constitutive_relation_name,
            std::unique_ptr<constitutive::BaseConstitutiveRelation>);
        void set_elements_constitutive_relation(ElementSelection selection,
            const std::string constitutive_relation_name);
        /// A boundary condition to the analysis
        /// Note: the BoundaryCondition object will store what constraints it creates based
        /// on the Node selection given to it when constructed.
        void add_boundary_condition(std::unique_ptr<BoundaryCondition> boundary_condition);

        /// Creates all of the objects needed before solve_time_steps can be called
        virtual void setup_analysis();
        /// Loops over the time steps, does iteration, writes output, etc.
        /// It is the primary analysis method.
        void solve_time_steps();
        /// Restart from a previously output restart file.
        void restart_from_file(const std::string filename);

    protected:
        std::unique_ptr<TimeStepManager> time_step_manager_;
        std::unique_ptr<Mesh> mesh_;
        flex_mesh::SubcomponentMaps subcomponent_maps_;
        std::unique_ptr<assembler::Assembler> assembler_;
        std::unique_ptr<dof_manager::DofMap> dof_map_;
        std::vector<std::unique_ptr<BoundaryCondition>> boundary_conditions_;
        std::map<std::string, std::unique_ptr<constitutive::BaseConstitutiveRelation>> constitutive_relations_;
		/// The global constraint matrix.  Right now, the interface for adding boundary conditions only allows
		/// for one constraint matrix.  So this property is put here instead of in a derived class.  If the
		/// interface is generalized to have multiple constraint matrices, then this will need to be moved.
		/// Note: this is still expected to be setup in initialize_linear_algebra in the derived classes.
		std::unique_ptr<sparse_linear_algebra::ConstraintMatrix> gl_constraint_matrix_;

        output::OutputFlagCache output_flag_cache_;
        std::vector<std::unique_ptr<output::OutputFlag>> output_flags_;
        fea_io::MetaWriter output_meta_;

        /// Initialize the global matrices and vectors for the analysis
		virtual void set_element_formulations() = 0;
        virtual void initialize_linear_algebra() = 0;
        virtual void apply_boundary_conditions();
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
            dof_map_->get_array_of_element_maps(),
            dof_map_->get_array_of_element_n_dofs(),
            static_cast<flex_mesh::LocalIndex>(mesh_->elements().size()),
            dof_map_->get_n_local_dofs(),
            dof_map_->get_n_global_dofs());
    }
}