#pragma once

#include "bvp/bc/BoundaryCondition.h"
#include "dof_manager/DofMap.h"
#include "finite_element/hhe_primitives/PrimitiveAlgorithms.h"
#include "mesh/common/Selection.h"
#include "sparse_linear_algebra/ConstraintMatrix.h"
#include "xfem/misc/GeneralDofUtilities.h"

namespace bvp
{
struct NeumannBCOpts
{
    /// Value to assign to the component of the field variable
    std::array<double, 3> value = { 0, 0, 0 };
    /// For the surface application type, the BC will be applied to the basal fields that touch the
    /// surface.  For the node application type, the BC will be applied DoFs cooresponding to the
    /// node where signs match up the enrichment tree (in other words, the node lies within the
    /// basal field's physical domain).
    BCApplicationType application_type = BCApplicationType::SURFACE;

    template<typename SerializerType>
    void serialize(SerializerType& s) const
    {
        s(utility::named_value("value", utility::slice(value)),
          utility::named_value("application_type", application_type));
    }

    template<typename SerializerType>
    void deserialize(SerializerType& s)
    {
        s(utility::slice(value), application_type);
    }
};

// TODO implement NeumannBC
}