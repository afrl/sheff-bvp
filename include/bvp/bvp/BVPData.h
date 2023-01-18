#pragma once

#include "assembler/Assembler.h"
#include "bvp/Defines.h"
#include "bvp/TimeDiscretization.h"
#include "bvp/bc/BoundaryCondition.h"
#include "bvp/bc/Serialization.h"
#include "dof_manager/DofManager.h"
#include "sparse_linear_algebra/SparseSolverOptions.h"
#include "xfem/XFEM.h"

namespace bvp
{
/// @brief Data that must be saved in a restart file
struct BVPCoreData
{
    std::string bvp_name;
    fs::path output_dir;
    TimeDiscretizationData time_discr_data;
    sparse_linear_algebra::SparseSolverOptions solver_opts;
    sparse_linear_algebra::NonLinearSolverOpts nonlinear_solver_opts;

    finite_element::FiniteElementMesh fe_mesh;
    dof_manager::DofManager dof_manager;
    std::vector<std::unique_ptr<BoundaryCondition>> boundary_conditions;
    /// Current solution vector from the last linear algebra solve for the owned and ghost DoFs
    /// on this rank. This is updated each solve, which means it may not be the converged solution
    /// for nonlinear problems. The order of DoFs follows the local indices of the DoF map.
    std::vector<double> solution;
    std::vector<xfem::EnrichmentSurfaceData> enrichment_surfaces_data;

    BVPCoreData() = default;
    BVPCoreData(BVPCoreData&&) noexcept;
};

/// @brief Data that can be recovered from BVPCoreData or set directly
struct BVPSecondaryData
{
    mesh::ElementAdjacencyList adjacency_list;
    mesh::NodeToElementMap node_to_element_map;
    dof_manager::ElementScatterMap element_scatter_map;
    std::unique_ptr<assembler::Assembler> assembler_ptr;
};

/// @brief Base class of struct that holds linear algebra objects (besides the solution vector),
/// allowing each type of BVP to define its own list of linear algebra structures.
struct BVPLinearAlgebra
{
    virtual ~BVPLinearAlgebra() {}

    virtual std::vector<sparse_linear_algebra::Vector*> vectors() const                       = 0;
    virtual std::vector<sparse_linear_algebra::SparseMatrix*> matrices() const                = 0;
    virtual std::vector<sparse_linear_algebra::ConstraintMatrix*> constraint_matrices() const = 0;
};

/// @brief Data structure that holds all of the data relevant to a boundary value problem.
///
/// For boundary value problems (BVP), the data (which this structure contains) is separated from
/// the algorithms.  All calls to a BoundaryValueProblem class require a BVPData instance, which
/// will be modified to solve the related boundary value problem.
struct BVPData
{
    /// Partially initialized as the BVP is set up, but fully initialized when
    /// BVP::solve_time_steps() is called
    BVPCoreData core_data;
    /// Initialized when BVP::solve_time_steps() is called
    BVPSecondaryData secondary_data;
    /// Initialized when BVP::solve_time_steps() is called
    std::unique_ptr<BVPLinearAlgebra> gl_lin_alg;

    BVPData() = default;
    BVPData(BVPData&&) noexcept;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Short-hand accessors to ease navigating the structs

    inline auto& bvp_name()
    {
        return core_data.bvp_name;
    }

    inline auto& bvp_name() const
    {
        return core_data.bvp_name;
    }

    inline auto& output_dir()
    {
        return core_data.output_dir;
    }

    inline auto& output_dir() const
    {
        return core_data.output_dir;
    }

    inline auto& mesh()
    {
        return core_data.fe_mesh.mesh;
    }

    inline auto& mesh() const
    {
        return core_data.fe_mesh.mesh;
    }

    inline auto& fe_mesh()
    {
        return core_data.fe_mesh;
    }

    inline auto& fe_mesh() const
    {
        return core_data.fe_mesh;
    }

    inline auto& dof_manager()
    {
        return core_data.dof_manager;
    }

    inline auto& dof_manager() const
    {
        return core_data.dof_manager;
    }

    inline auto& bcs()
    {
        return core_data.boundary_conditions;
    }

    inline auto& solution()
    {
        _UTIL_ASSERT_C(core_data.solution.size() > 0,
                       "Linear algebra must be intiailized before the solution can be accessed.");
        return core_data.solution;
    }

    inline auto& solution() const
    {
        _UTIL_ASSERT_C(core_data.solution.size() > 0,
                       "Linear algebra must be intiailized before the solution can be accessed.");
        return core_data.solution;
    }

    inline auto& enrichment_surfaces_data()
    {
        return core_data.enrichment_surfaces_data;
    }

    inline auto& adjacency_list()
    {
        return secondary_data.adjacency_list;
    }

    inline auto& node_to_element_map()
    {
        return secondary_data.node_to_element_map;
    }

    inline auto& element_scatter_map()
    {
        return secondary_data.element_scatter_map;
    }

    inline auto& assembler()
    {
        _UTIL_ASSERT_C(secondary_data.assembler_ptr,
                       "The assembler must be initialized first, BoundaryValueProblem is "
                       "probably missing a call to init_data internally.");
        return *secondary_data.assembler_ptr;
    }

    inline void add_bc(std::unique_ptr<BoundaryCondition>& bc)
    {
        core_data.boundary_conditions.push_back(std::move(bc));
    }

    template<typename DerivedBC>
    inline void add_bc(DerivedBC&& bc)
    {
        core_data.boundary_conditions.push_back(std::make_unique<DerivedBC>(std::move(bc)));
    }
};

/// @brief Queries how many ranks were running for BVPData saved to the given file.
/// @param filestem - name of the file to load (extension will be ignored if given)
/// @return Number of ranks that saved the data.
int query_n_ranks_for_bvp_data(const std::string filestem);

/// @brief Load a binary serialization of a BVPData object in the given file.
/// @param filestem - name of the file to load (extension will be ignored if given)
///
/// Note: if the BVPData was saved while distributed, then it is required that the current number of
/// ranks match how many were running when it was saved.  In this case, this loads the data
/// corresponding to the current rank ID.
BVPData load_bvp_data(const std::string filestem);

/// @brief Load a binary serialization of a BVPData object in the given file.
/// @param filestem - name of the file to load (extension will be ignored if given)
/// @param rank_to_load - the corresponding rank for which to load the saved BVPData
///
/// Note: As opposed to the other load_bvp_data calls, this does not require the current number of
/// ranks to match the number that saved the data and the current rank does not have to match the
/// rank for which data will be loaded.
BVPData load_bvp_data(const std::string filestem, int rank_to_load);

/// @brief Write a binary serialization of the given bvp_data to a file with the given name and a
/// corresponding JSON file containing the offsets of each rank.
/// @param bvp_data - the BVPData to save
/// @param filestem - name of the file without an extension (can be a path relative to the CWD)
///
/// Note: writes the binary file to "<filestem>.bvp.bin" and the JSON data to "<filestem>.bvp.json"
void save_bvp_data(const BVPData& bvp_data, const std::string filestem);

/// @brief Write a JSON serialization of the given bvp_data to a text file with the given name.
/// @param bvp_data - data to save
/// @param filestem - name of the file without an extension (can be a path relative to the CWD)
///
/// Note: if distributed, then each rank will write to its own text file "<filestem>_<rank id>.txt",
/// otherwise this will write to  "<filestem>.txt"
void save_bvp_data_json(const BVPData& bvp_data, const std::string filestem);

/// @brief Structures in the sparse_linear_algebra modules require references to DoF maps, and this
/// function creates the structure holding the references from the dof_manager.  It is called by the
/// BVPCoreData move constructor and BoundaryValueProblem::init_core_data.
sparse_linear_algebra::DofMapRefs create_dof_refs(const dof_manager::DofManager& dof_manager);
}
