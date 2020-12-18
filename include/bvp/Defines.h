#pragma once

#include "finite_element/FiniteElement.h"
#include "mesh/Mesh.h"

namespace bvp
{
/// Must be called ONLY once and BEFORE an analysis is reaches the setup phase.
/// \param preferred_lin_alg_backend - a string referring to the type of backend
///        to use for linear algebra. (Ex: MKL, PETSC, PARLAB)
/// \param argc - count of commandline arguments to pass on to other modules.
/// \param argc - commandline arguments to pass on to other modules.
void init(const std::string preferred_lin_alg_backend, int* argc, char*** argv);

/// @brief A simple struct to hold global options for the BVP module.
/// Typically these should not be modified except for debugging and unit tests.
struct BVPGlobalOptions
{
    static BVPGlobalOptions& get()
    {
        static BVPGlobalOptions opts;
        return opts;
    }

    bool save_matrices = false;

  private:
    BVPGlobalOptions() = default;
};
}