#pragma once

#include "bvp/bc/BoundaryCondition.h"
#include "dof_manager/DofMap.h"
#include "finite_element/hhe_primitives/PrimitiveAlgorithms.h"
#include "mesh/common/Selection.h"
#include "xfem/misc/GeneralDofUtilities.h"

namespace bvp
{
struct TractionBCOpts
{
    /// Value to assign to the component of the field variable
    simple_math::Vector3D traction = { 0, 0, 0 };

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("traction", utility::slice(traction.data(), traction.size())));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(utility::slice(traction.data(), traction.size()));
    }
};

class TractionBC final : public BoundaryCondition
{
  public:
    TractionBC() {}
    TractionBC(const std::string& node_group_name, const TractionBCOpts& opts)
      : node_group_name_(node_group_name)
      , opts_(opts)
    {
    }
    virtual ~TractionBC() {}

    /// @see BoundaryCondition::type
    virtual BCType type() const override
    {
        return BCType::TRACTION;
    }

    /// @see BoundaryCondition::update_constraints
    virtual void update_constraints() override {}

    /// @see BoundaryCondition::initialize
    virtual void initialize(const finite_element::FiniteElementMesh& fe_mesh,
                            const mesh::NodeToElementMap& node_to_element_map,
                            const dof_manager::DofMap& dof_map,
                            const dof_manager::RankDofIndex& n_components_for_field,
                            const double& scale_factor,
                            const bool& is_incremental,
                            sparse_linear_algebra::ConstraintMatrix& constraint_matrix) override;

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_ref("node_group_name", node_group_name_),
          utility::named_ref("opts", opts_));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(node_group_name_, opts_);
    }

  protected:
    std::string node_group_name_;
    TractionBCOpts opts_;
};
}