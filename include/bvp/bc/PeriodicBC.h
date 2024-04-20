#pragma once

#include "bvp/bc/BoundaryCondition.h"

namespace bvp
{

class PeriodicBC final : public BoundaryCondition
{
  public:
    PeriodicBC() {}
    PeriodicBC(const std::string& primary_node_group_name,
               const std::string& secondary_node_group_name,
               const mesh::Point& periodicity_vector,
               const BCApplicationType& application_type = BCApplicationType::SURFACE)
      : primary_node_group_name_(primary_node_group_name)
      , secondary_node_group_name_(secondary_node_group_name)
      , periodicity_vector_(periodicity_vector)
      , application_type_(application_type)
    {
    }
    virtual ~PeriodicBC() {}

    /// @see BoundaryCondition::type
    virtual BCType type() const override
    {
        return BCType::PERIODIC;
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

    /// @see BoundaryCondition::get_scale_factor_function_index
    virtual int get_scale_factor_function_index() const override
    {
        return -1;
    }

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_ref("primary_node_group_name", primary_node_group_name_),
          utility::named_ref("secondary_node_group_name", secondary_node_group_name_),
          utility::named_ref("periodicity_vector", periodicity_vector_),
          utility::named_value("application_type", application_type_));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(primary_node_group_name_,
          secondary_node_group_name_,
          periodicity_vector_,
          application_type_);
    }

  protected:
    /// Nodes that correspond to the primary DoFs in the multipoint constraints
    std::string primary_node_group_name_;
    /// Nodes that correspond to the secondary DoFs, which are a function of the primary DoFs
    std::string secondary_node_group_name_;
    /// Value to assign to the component of the field variable, default 0
    mesh::Point periodicity_vector_ = { 0., 0., 0. };
    /// For the surface application type, the BC will be applied to the basal fields that touch the
    /// surface.  For the node application type, the BC will be applied DoFs cooresponding to the
    /// node where signs match up the enrichment tree (in other words, the node lies within the
    /// basal field's physical domain).
    BCApplicationType application_type_ = BCApplicationType::SURFACE;
};
}