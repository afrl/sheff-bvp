#include "bvp/bc/TractionBC.h"

#include "finite_element/calc/IntegrateOverFace.h"
#include "mesh/algorithms/Search.h"

namespace bvp
{
void TractionBC::initialize(const finite_element::FiniteElementMesh& fe_mesh,
                            const mesh::NodeToElementMap& node_to_element_map,
                            const dof_manager::DofMap& dof_map,
                            const dof_manager::RankDofIndex& n_components_for_field,
                            const double& scale_factor,
                            const bool& is_incremental,
                            sparse_linear_algebra::ConstraintMatrix& constraint_matrix)
{
    std::set<dof_manager::RankDofIndex> dofs_to_constrain;
    const auto& rank_to_global_map = dof_map.rank_to_global_map;

    std::set<mesh::LocalIndex> node_indices_set(
      fe_mesh.mesh.nodes().group(node_group_name_).node_indices.begin(),
      fe_mesh.mesh.nodes().group(node_group_name_).node_indices.end());
    auto selected_faces = mesh::find_faces_from_node_set(fe_mesh.mesh, node_indices_set);

    std::vector<std::pair<sparse_linear_algebra::LocalIndex, double>> dof_force_terms_to_add;
    for (auto& face : selected_faces)
    {
        auto fe_accessor    = fe_mesh.get_finite_element(face.local_element_index);
        auto n_basal_fields = fe_accessor.layout.get_n_basal_fields();
        xfem::ElementXFEMDofMap element_dof_map(dof_map, n_basal_fields, fe_accessor);

        auto vector_function = [traction = opts_.traction](const simple_math::Vector3D&) {
            return traction;
        };
        auto integration_result =
          finite_element::calculate_integrated_nodal_values_for_vector_function(fe_accessor,
                                                                                face.face_index,
                                                                                0,
                                                                                vector_function);

        for (auto& patch_integration_result : integration_result)
        {
            auto& bounded_basal_field_index = patch_integration_result.first;
            auto& integrated_nodal_values   = patch_integration_result.second;
            auto bounded_field_dof_map      = element_dof_map.at_leaf(bounded_basal_field_index);

            for (auto i : utility::range(fe_accessor.element.connectivity().size()))
            {
                for (auto j : utility::range(3))
                {
                    if (integrated_nodal_values(i, j) != 0)
                    {
                        dof_force_terms_to_add.push_back(
                          { bounded_field_dof_map.at_node(i)[j], integrated_nodal_values(i, j) });
                    }
                }
            }
        }
    }

    for (auto& dof_force_term_to_add : dof_force_terms_to_add)
    {

        double force_value = scale_factor * dof_force_term_to_add.second;

        const auto& global_dof = rank_to_global_map[dof_force_term_to_add.first];
        constraint_matrix.add_neumann_term(global_dof, force_value);
    }
}
}