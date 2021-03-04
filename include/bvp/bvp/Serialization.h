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
    s(utility::named_value("output_dir", bvp_data.core_data.output_dir.generic_string()),
      utility::named_ref("time_discr_data", bvp_data.core_data.time_discr_data),
      utility::named_ref("fe_mesh", bvp_data.core_data.fe_mesh),
      utility::named_ref("dof_manager", bvp_data.core_data.dof_manager),
      utility::named_value("n_boundary_conditions", bvp_data.core_data.boundary_conditions.size()),
      utility::named_value("boundary_conditions",
                           utility::slice(bvp_data.core_data.boundary_conditions.data(),
                                          bvp_data.core_data.boundary_conditions.size())));

    s(utility::named_value("has_gl_solution", bvp_data.core_data.gl_solution != nullptr));
    if (bvp_data.core_data.gl_solution)
    {
        s(utility::named_ref("gl_solution", *bvp_data.core_data.gl_solution));
    }

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
    s(output_dir_str,
      bvp_data.core_data.time_discr_data,
      bvp_data.core_data.fe_mesh,
      bvp_data.core_data.dof_manager);
    bvp_data.core_data.output_dir = output_dir_str;

    std::size_t n_boundary_conditions;
    s(n_boundary_conditions);
    bvp_data.core_data.boundary_conditions.resize(n_boundary_conditions);
    s(utility::slice(bvp_data.core_data.boundary_conditions.data(),
                     bvp_data.core_data.boundary_conditions.size()));

    bool has_gl_solution;
    s(has_gl_solution);
    if (has_gl_solution)
    {
        auto& dof_maps                             = bvp_data.dof_manager().get_dof_maps();
        sparse_linear_algebra::DofMapRefs dof_refs = {
            dof_maps.n_owned_dofs,           dof_maps.n_local_ghost_dofs,
            dof_maps.n_exclusive_ghost_dofs, dof_maps.owned_global_dof_begin,
            dof_maps.owned_global_dof_end,   dof_maps.n_owned_elements,
            dof_maps.element_to_rank_maps,   dof_maps.element_n_dofs,
            dof_maps.rank_to_global_map
        };

        bvp_data.core_data.gl_solution =
          sparse_linear_algebra::LinAlgFactory::create_vector(dof_refs);
        bvp_data.core_data.gl_solution->allocate();

        s(*bvp_data.core_data.gl_solution);
    }

    std::size_t n_enrichment_surfaces = 0;
    s(n_enrichment_surfaces);
    bvp_data.core_data.enrichment_surfaces_data.resize(n_enrichment_surfaces);
    s(utility::slice(bvp_data.core_data.enrichment_surfaces_data.data(),
                     bvp_data.core_data.enrichment_surfaces_data.size()));
}
}