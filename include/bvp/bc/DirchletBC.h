#pragma once

#include "bvp/bc/BoundaryCondition.h"

namespace bvp
{
struct DirichletBCOpts
{
    /// Component of the field variable, default 0
    unsigned short field_component = 0;
    /// Value to assign to the component of the field variable, default 0
    double value = 0;
    /// For the surface application type, the BC will be applied to the basal fields that touch the
    /// surface.  For the node application type, the BC will be applied DoFs cooresponding to the
    /// node where signs match up the enrichment tree (in other words, the node lies within the
    /// basal field's physical domain).
    BCApplicationType application_type = BCApplicationType::SURFACE;

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_ref("field_component", field_component),
          utility::named_ref("value", value),
          utility::named_value("application_type", application_type));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(field_component, value, application_type);
    }
};

class DirichletBC final : public BoundaryCondition
{
  public:
    DirichletBC() {}
    DirichletBC(const std::string& node_group_name, const DirichletBCOpts& opts)
      : node_group_name_(node_group_name)
      , opts_(opts)
    {
    }
    virtual ~DirichletBC() {}

    /// @see BoundaryCondition::type
    virtual BCType type() const override
    {
        return BCType::DIRICHLET;
    }

    /// @see BoundaryCondition::update_constraints
    virtual void update_constraints() override {}

    /// @see BoundaryCondition::initialize
    virtual void initialize(const finite_element::FiniteElementMesh& fe_mesh,
                            const mesh::NodeToElementMap& node_to_element_map,
                            const dof_manager::DofMap& dof_map,
                            const dof_manager::RankDofIndex& n_components_for_field,
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
    DirichletBCOpts opts_;
};
}