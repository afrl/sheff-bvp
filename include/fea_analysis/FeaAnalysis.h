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

        /// @name Pre-Solve Methods
        /// Methods that should be called before solve_time_steps() is called.
        ///@{

        /// Sets the input for the analysis, reads the input mesh (and any groups or data
        /// assocciated with it), and sets up the writer for the output meta file.
        /// \param input_script_path - the path to input script is used for information
        ///        in the output meta file.  The file is not accessed by this method.
        /// \param input_meta_file_path - the path to the input meta file that contains the
        ///        the mesh, groups, and input data.  For an analysis from the start, this
        ///        should be the path to the file created by the preprocessor, but for a
        ///        restart, it should be the output meta file that is being restarted.
        ///        Note: the mesh, groups, and input data will be copied to the output
        ///        directory.  This ensures that the output directory is self-contained, but
        ///        does come at a cost of hard disk space.
        /// \param output_dir_path - the path to the output directory.  This can be a relative
        ///        path but should be relative to the current working directory.  If the
        ///        directory does not exist yet, it will be created.
        void set_input(
            const std::string input_script_path,
            const std::string input_meta_file_path,
            const std::string output_dir_path);
        /// Sets input for the analysis, similar to the other set_input() method, but
        /// this method sets the mesh object directory and is intended for unit tests
        /// or scripts that create the mesh on the fly.
        /// \param mesh - the input mesh to use for the analysis
        /// \param output_dir_path - the path to the output directory.  This can be a relative
        ///        path but should be relative to the current working directory.  If the
        ///        directory does not exist yet, it will be created.
        void set_input(
            std::unique_ptr<Mesh> mesh,
            const std::string analysis_name,
            const std::string output_dir_path);
        /// Add a constitutive relation to the analysis
        void add_constitutive_relation(const std::string constitutive_relation_name,
            std::unique_ptr<constitutive::BaseConstitutiveRelation>);
        void set_elements_constitutive_relation(ElementSelection selection,
            const std::string constitutive_relation_name);
        /// A boundary condition to the analysis
        /// Note: the BoundaryCondition object will store what constraints it creates based
        /// on the Node selection given to it when constructed.
        void add_boundary_condition(std::unique_ptr<BoundaryCondition> boundary_condition);

        /// Creates all of the objects needed before solve_time_steps() can be called.
        /// This should be called after all other pre-solve methods have been called, but
        /// before pre-cracks are inserted (if applicable).
        virtual void setup_analysis();

        ///@}
        // End Pre-Solve Methods

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
        std::unique_ptr<sparse_linear_algebra::DofMapReferences> dof_refs_;

        output::OutputFlagCache output_flag_cache_;
        std::vector<std::unique_ptr<output::OutputFlag<finite_element::VolumeQuantity>>> volume_output_flags_;
        std::vector<std::unique_ptr<output::OutputFlag<finite_element::SurfaceQuantity>>> surface_output_flags_;
        fea_io::MetaWriter output_meta_;
        fs::path output_dir_;
        fs::path curr_time_step_dir_;

		virtual void set_element_formulations() = 0;
        virtual void initialize_linear_algebra() = 0;
        virtual void apply_boundary_conditions();
        virtual void solve_time_step() = 0;
        virtual void write_output() = 0;
        virtual void write_restart_file() = 0;
    };
}