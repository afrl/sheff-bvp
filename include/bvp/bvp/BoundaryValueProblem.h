#pragma once

#include "BVPData.h"
#include "bvp/Defines.h"
#include "finite_element/hhe_primitives/BasalFieldSelection.h"

#include <vector>

namespace bvp
{

struct SubstepResult
{
    bool converged = false;
    double time    = 0.0;
    sparse_linear_algebra::SolverStatus solver_status;
};

// @brief Simple class describing a result of calling BoundaryValueProblem::solve_next_time_step.
struct SolveTimeStepResult
{
    bool converged         = false;
    bool end_of_time_steps = false;
    std::vector<SubstepResult> substep_results;
};

/// @brief General base class for all types of boundary value problems.
///
/// Follows singleton pattern, so all derived types should implement a method like "static
/// DerivedBVP& get()".
class BoundaryValueProblem
{
  public:
    virtual ~BoundaryValueProblem() {}

    /// @name Pre-Solve Methods
    /// Methods that should be called before solve_time_steps() is called.
    ///@{

    /// Sets the input for the analysis, reads the input mesh (and any groups or data
    /// assocciated with it).
    /// @param input_script_path - the path to input script is used for information
    ///        in the output meta file.  The file is not accessed by this method.
    /// @param input_meta_file_path - the path to the input meta file that contains the
    ///        the mesh, groups, and input data.  For an analysis from the start, this
    ///        should be the path to the file created by the preprocessor, but for a
    ///        restart, it should be the output meta file that is being restarted.
    ///        Note: the mesh, groups, and input data will be copied to the output
    ///        directory.  This ensures that the output directory is self-contained, but
    ///        does come at a cost of hard disk space.
    /// @param output_dir_path - the path to the output directory.  This can be a relative
    ///        path but should be relative to the current working directory.  If the
    ///        directory does not exist yet, it will be created.
    /// @param bvp_data - BVP data to modify
    void set_input(const std::string input_script_path,
                   const std::string input_meta_file_path,
                   const std::string output_dir_path,
                   BVPData& bvp_data) const;

    /// Sets input for the analysis, similar to the other set_input() method, but
    /// this method sets the mesh object directory and is intended for unit tests
    /// or scripts that create the mesh on the fly.
    /// @param mesh - the input mesh to use for the analysis
    /// @param analysis_name - a name for the analysis (used for post-processing)
    /// @param output_dir_path - the path to the output directory.  This can be a relative
    ///        path but should be relative to the current working directory.  If the
    ///        directory does not exist yet, it will be created.
    /// @param bvp_data - BVP data to modify
    void set_input(mesh::Mesh&& mesh,
                   const std::string analysis_name,
                   const std::string output_dir_path,
                   BVPData& bvp_data) const;

    /// @brief Sets the model for basal fields based on the formulation type and properties.
    ///
    /// Setting the model includes the point-wise model, the finite element model, and
    /// all model data.  The combination of objects that compose the model is determined
    /// based on the combination fo the formulation and properties given.
    ///
    /// @param selection - the selection object that will return the basal fields to assign
    ///        a field model to based on the other given options.
    /// @param formulation_name - valid options for formulations currently include:
    ///        - "InfinitesimalStrain"
    /// @param properties_name - the properties to be assigned to the basal fields.  This
    ///        analysis object makes a deep copy of the properties to ensure it will
    ///        be persistent throughout the analysis.  A copy is created for each call,
    ///        so be cognizant that calling this method many times with the same properties
    ///        object will be far less memory-efficient than calling it once with a more
    ///        complex selection.
    /// @param bvp_data - BVP data to modify
    void set_field_models(finite_element::BasalFieldSelection selection,
                          const std::string& formulation_name,
                          const std::string& properties_name,
                          BVPData& bvp_data) const;

    /// @brief Sets the model for basal fields based on the formulation type and properties.
    ///
    /// See other signature for more documentation.  This provides the same functionality except
    /// that the list of basal fields are explicitly given.
    void set_field_models(const std::vector<finite_element::SelectedBasalField>& selected_fields,
                          const std::string& formulation_name,
                          const std::string& properties_name,
                          BVPData& bvp_data) const;

    ///@}
    // End Pre-Solve Methods

    /// @brief Solves all of the time steps for the boundary value problem.
    ///
    /// Note: equivalent to calling solve_next_time_step until there are no more time steps.
    void solve_time_steps(BVPData& bvp_data) const;

    /// @brief Solves the next time step for the boundary value problem, if there any remaining.
    virtual SolveTimeStepResult solve_next_time_step(BVPData& bvp_data) const = 0;

  protected:
    BoundaryValueProblem() {}

    /// @brief Ensures the secondary data and linear algebra is initialized, which is required
    /// before solves take place or enrichment surfaces are inserted.
    ///
    /// Required to be called before a solve takes place; however, it is automatically called within
    /// other methods of this class (and derived classes) when needed.
    virtual void init_all_data(BVPData& bvp_data) const;

    /// @brief Ensures all of the core data is initialized.
    ///
    /// Required to be called before a solve takes place; however, it is automatically called within
    /// other methods of this class (and derived classes) when needed.
    void init_core_data(BVPData& bvp_data) const;

    /// @brief Initializes the linear algebra containers (bvp_data.core_data.gl_solution and
    /// bvp_data.gl_lin_alg).
    ///
    /// Sets bvp_data.gl_lin_alg.initialized to true.
    virtual void init_lin_algebra(BVPData& bvp_data) const = 0;

    /// @brief Applies the boundary conditions in bvp_data.boundary_conditions to the system
    /// of equations
    ///
    /// Modifies the constraint matrix and appropriate vectors in bvp_data.gl_lin_alg.
    virtual void apply_boundary_conditions(BVPData& bvp_data) const = 0;
};
}