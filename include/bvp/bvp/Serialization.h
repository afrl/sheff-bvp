#pragma once

#include "BVPData.h"
#include "sparse_linear_algebra/LinAlgFactory.h"
#include "utility/Slice.h"
#include "utility/serialization/NamedValue.h"
#include "xfem/XFEM.h"

namespace utility
{
/// @brief Provides serialization support for BVPData
template<typename SerializerType>
void serialize(SerializerType& s, const bvp::BVPData& bvp_data)
{
    s(utility::named_value("bvp_name", bvp_data.core_data.bvp_name),
      utility::named_value("output_dir", bvp_data.core_data.output_dir.generic_string()),
      utility::named_ref("time_discr_data", bvp_data.core_data.time_discr_data),
      utility::named_ref("solver_opts", bvp_data.core_data.solver_opts),
      utility::named_ref("fe_mesh", bvp_data.core_data.fe_mesh),
      utility::named_ref("dof_manager", bvp_data.core_data.dof_manager),
      utility::named_value("n_boundary_conditions", bvp_data.core_data.boundary_conditions.size()),
      utility::named_value("boundary_conditions",
                           utility::slice(bvp_data.core_data.boundary_conditions.data(),
                                          bvp_data.core_data.boundary_conditions.size())));
    // Since we are not writing the solution size, assert it is the expected size to avoid seg
    // faults if something goes wrong.
    assert(bvp_data.core_data.dof_manager.get_dof_maps().n_dofs_on_rank() ==
           bvp_data.core_data.solution.size());
    s(utility::named_value("solution", utility::slice(bvp_data.core_data.solution)));
    s(utility::named_value("n_enrichment_surfaces",
                           bvp_data.core_data.enrichment_surfaces_data.size()),
      utility::named_value("enrichment_surfaces_data",
                           utility::slice(bvp_data.core_data.enrichment_surfaces_data.data(),
                                          bvp_data.core_data.enrichment_surfaces_data.size())));
}

/// @brief Provides deserialization support for boundary conditions (since they are polymorphic).
///
/// Reads the BC type, creates an a new object of the appropriate type, and forwards the call to the
/// concrete class's deserialize method.
template<typename SerializerType>
void deserialize(SerializerType& s, bvp::BVPData& bvp_data)
{
    std::string output_dir_str;
    s(bvp_data.core_data.bvp_name,
      output_dir_str,
      bvp_data.core_data.time_discr_data,
      bvp_data.core_data.solver_opts,
      bvp_data.core_data.fe_mesh,
      bvp_data.core_data.dof_manager);
    bvp_data.core_data.output_dir = output_dir_str;

    std::size_t n_boundary_conditions;
    s(n_boundary_conditions);
    bvp_data.core_data.boundary_conditions.resize(n_boundary_conditions);
    s(utility::slice(bvp_data.core_data.boundary_conditions.data(),
                     bvp_data.core_data.boundary_conditions.size()));

    bvp_data.core_data.solution.resize(bvp_data.dof_manager().get_dof_maps().n_dofs_on_rank());
    s(utility::slice(bvp_data.core_data.solution));

    std::size_t n_enrichment_surfaces = 0;
    s(n_enrichment_surfaces);
    bvp_data.core_data.enrichment_surfaces_data.resize(n_enrichment_surfaces);
    s(utility::slice(bvp_data.core_data.enrichment_surfaces_data.data(),
                     bvp_data.core_data.enrichment_surfaces_data.size()));
}
}