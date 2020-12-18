#pragma once

#include "BoundaryCondition.h"
#include "DirchletBC.h"
#include "utility/serialization/NamedValue.h"

#include <memory>

namespace utility
{
/// @brief Provides serialization support for boundary conditions (since they are polymorphic).
///
/// Writes the BC type and forwards the call to the concrete class's serialize method.
template<typename SerializerType>
void serialize(SerializerType& s, const std::unique_ptr<bvp::BoundaryCondition>& bc)
{
    bvp::BCType type = bc->type();
    s(utility::named_value("type", type));

    if (type == bvp::BCType::DIRICHLET)
    {
        auto& typed_bc = static_cast<const bvp::DirichletBC&>(*bc);
        s(utility::named_ref("data", typed_bc));
        return;
    }
    // TODO implement case for NeumannBC once it is implemented
    else
    {
        throw _UTIL_MSG_EXCEPT("Case not implemented to forward serialization for given BC type: " +
                               bvp::to_string(type));
    }
}

/// @brief Provides deserialization support for boundary conditions (since they are polymorphic).
///
/// Reads the BC type, creates an a new object of the appropriate type, and forwards the call to the
/// concrete class's deserialize method.
template<typename SerializerType>
void deserialize(SerializerType& s, std::unique_ptr<bvp::BoundaryCondition>& bc)
{
    bvp::BCType type;
    s(type);

    if (type == bvp::BCType::DIRICHLET)
    {
        bc             = std::make_unique<bvp::DirichletBC>();
        auto& typed_bc = static_cast<bvp::DirichletBC&>(*bc);
        s(typed_bc);
        return;
    }
    // TODO implement case for NeumannBC once it is implemented
    else
    {
        throw _UTIL_MSG_EXCEPT(
          "Case not implemented to forward deserialization for given BC type: " +
          bvp::to_string(type));
    }
}
}