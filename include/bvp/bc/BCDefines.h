#pragma once

#include "utility/Exceptions.h"

#include <cinttypes>
#include <string>

namespace bvp
{
/// @brief Enumeration of all concrete derived types of BoundaryCondition.
enum class BCType : int8_t
{
    DIRICHLET,
    NEUMANN,
    TRACTION,
    PERIODIC
};

static std::string to_string(const BCType& type)
{
    switch (type)
    {
        case BCType::DIRICHLET: return "DIRICHLET";
        case BCType::NEUMANN: return "NEUMANN";
        case BCType::TRACTION: return "TRACTION";
        case BCType::PERIODIC: return "PERIODIC";
        default: throw _UTIL_MSG_EXCEPT("String conversion not implemented for given BCType enum.");
    }
}

/// @brief Enumeration of types of boundary applications that are supported.
enum class BCApplicationType : int8_t
{
    NODE,   ///< Boundary conditions are applied to specified nodes
    SURFACE ///< Boundary conditions are only applied faces that compose the surface
};

static std::string to_string(const BCApplicationType& type)
{
    switch (type)
    {
        case BCApplicationType::NODE: return "NODE";
        case BCApplicationType::SURFACE: return "SURFACE";
        default:
            throw _UTIL_MSG_EXCEPT(
              "String conversion not implemented for given BCApplicationType enum.");
    }
}
}