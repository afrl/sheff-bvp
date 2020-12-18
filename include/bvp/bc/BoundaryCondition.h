#pragma once

#include "BCDefines.h"
#include "dof_manager/DofMap.h"
#include "finite_element/FiniteElementMesh.h"
#include "sparse_linear_algebra/ConstraintMatrix.h"

namespace bvp
{
/// @brief Defines the interface for all boundary condition types.
class BoundaryCondition
{
  public:
    BoundaryCondition() {}
    ~BoundaryCondition() {}

    /// @brief Returns an enumeraiton indicating what concrete type of BC this object is.
    ///
    /// Required for serialization and deserialization.
    virtual BCType type() const = 0;
    /// TODO: assess if this method is actually needed.
    virtual void update_constraints() = 0;
    /// Called by BoundaryValueProblem::apply_boundary_conditions
    virtual void initialize(const finite_element::FiniteElementMesh& fe_mesh,
                            const mesh::NodeToElementMap& node_to_element_map,
                            const dof_manager::DofMap& dof_map,
                            const dof_manager::RankDofIndex& n_components_for_field,
                            sparse_linear_algebra::ConstraintMatrix& constraint_matrix) = 0;
};
}