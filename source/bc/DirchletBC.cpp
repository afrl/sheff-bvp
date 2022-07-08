#include "bvp/bc/DirchletBC.h"

#include "finite_element/hhe_primitives/PrimitiveAlgorithms.h"
#include "mesh/algorithms/Search.h"
#include "xfem/misc/GeneralDofUtilities.h"

namespace bvp
{
void DirichletBC::initialize(const finite_element::FiniteElementMesh& fe_mesh,
                             const mesh::NodeToElementMap& node_to_element_map,
                             const dof_manager::DofMap& dof_map,
                             const dof_manager::RankDofIndex& n_components_for_field,
                             const double& scale_factor,
                             const bool& is_incremental,
                             sparse_linear_algebra::ConstraintMatrix& constraint_matrix)
{
    std::set<dof_manager::RankDofIndex> dofs_to_constrain;
    const auto& rank_to_global_map = dof_map.rank_to_global_map;

    if (opts_.application_type == BCApplicationType::NODE)
    {
        auto& nodes_constrain = fe_mesh.mesh.nodes().group(node_group_name_).node_indices;
        for (auto& node_index : nodes_constrain)
        {
            for (auto& parent_element : node_to_element_map.get_elements(node_index))
            {
                auto fe_accessor    = fe_mesh.get_finite_element(parent_element.local_index());
                auto n_basal_fields = fe_accessor.layout.get_n_basal_fields();
                xfem::ElementXFEMDofMap element_dof_map(dof_map, n_basal_fields, fe_accessor);
                finite_element::EnrichmentTreeTraverser traverser(fe_accessor);

                mesh::ConnIndex node_conn_index = std::numeric_limits<mesh::ConnIndex>::max();
                for (auto i : utility::range(parent_element.connectivity().size()))
                {
                    if (parent_element.connectivity()[i] == node_index)
                    {
                        node_conn_index = i;
                    }
                }
                assert(node_conn_index != std::numeric_limits<mesh::ConnIndex>::max());

                for (auto basal_field_index : utility::range(n_basal_fields))
                {
                    auto& basal_field = fe_accessor.layout.get_basal_field(basal_field_index);
                    auto signs_match_result =
                      traverser.get_node_field_sign_equivalence_vector(basal_field.field_index);

                    if (signs_match_result[node_conn_index])
                    {
                        auto field_dof_map = element_dof_map.at_leaf(basal_field_index);
                        dofs_to_constrain.insert(
                          field_dof_map.at_node(node_conn_index)[opts_.field_component]);
                    }
                }
            }
        }
    }
    else if (opts_.application_type == BCApplicationType::SURFACE)
    {
        std::set<mesh::LocalIndex> node_indices_set(
          fe_mesh.mesh.nodes().group(node_group_name_).node_indices.begin(),
          fe_mesh.mesh.nodes().group(node_group_name_).node_indices.end());
        auto selected_faces = mesh::find_faces_from_node_set(fe_mesh.mesh, node_indices_set);

        for (auto& face : selected_faces)
        {
            auto fe_accessor    = fe_mesh.get_finite_element(face.local_element_index);
            auto n_basal_fields = fe_accessor.layout.get_n_basal_fields();
            xfem::ElementXFEMDofMap element_dof_map(dof_map, n_basal_fields, fe_accessor);
            mesh::FaceContainer faces(fe_accessor.element);

            for (auto basal_field_index : utility::range(n_basal_fields))
            {
                auto field_touches_face = finite_element::does_basal_field_touch_topology(
                  fe_accessor,
                  basal_field_index,
                  { mesh::TopologyType::FACE, face.face_index });
                if (field_touches_face)
                {
                    auto field_dof_map = element_dof_map.at_leaf(basal_field_index);

                    // Loop over nodes in face and add rank DoF's to set
                    for (auto i : faces[face.face_index].element_conn_indices())
                    {
                        dofs_to_constrain.insert(field_dof_map.at_node(i)[opts_.field_component]);
                    }
                }
            }
        }
    }

    double constraint_value = scale_factor * opts_.value;

    for (auto rank_dof : dofs_to_constrain)
    {
        const auto& global_dof = rank_to_global_map[rank_dof];
        constraint_matrix.add_dirichlet_constraint(global_dof, constraint_value);
    }
}
}