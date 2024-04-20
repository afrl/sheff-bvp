#include "bvp/bc/PeriodicBC.h"

#include "finite_element/hhe_primitives/PrimitiveAlgorithms.h"
#include "mesh/algorithms/Search.h"
#include "mesh/algorithms/Sort.h"
#include "xfem/misc/GeneralDofUtilities.h"

namespace bvp
{
static std::vector<std::vector<dof_manager::RankDofIndex>>
find_dofs_to_constrain_by_field_component(const std::string& node_group_name,
                                          const BCApplicationType& application_type,
                                          const std::vector<unsigned short> field_components,
                                          const finite_element::FiniteElementMesh& fe_mesh,
                                          const mesh::NodeToElementMap& node_to_element_map,
                                          const dof_manager::DofMap& dof_map)
{
    std::vector<std::vector<dof_manager::RankDofIndex>> dofs_to_constrain_by_field_component;
    dofs_to_constrain_by_field_component.resize(field_components.size());
    if (application_type == BCApplicationType::NODE)
    {
        auto nodes_to_constrain = fe_mesh.mesh.nodes().group(node_group_name).node_indices;
        mesh::sort_node_indices_by_coord(fe_mesh.mesh, nodes_to_constrain);
        for (auto& node_index : nodes_to_constrain)
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
                        for (auto i : utility::indices(field_components))
                        {
                            auto& field_component = field_components[i];
                            dofs_to_constrain_by_field_component[i].push_back(
                              field_dof_map.at_node(node_conn_index)[field_component]);
                        }
                    }
                }
            }
        }
    }
    else if (application_type == BCApplicationType::SURFACE)
    {
        std::set<mesh::LocalIndex> node_indices_set(
          fe_mesh.mesh.nodes().group(node_group_name).node_indices.begin(),
          fe_mesh.mesh.nodes().group(node_group_name).node_indices.end());
        auto selected_faces = mesh::find_faces_from_node_set(fe_mesh.mesh, node_indices_set);
        mesh::sort_element_face_pairs_by_coord(fe_mesh.mesh, selected_faces);

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
                        for (auto j : utility::indices(field_components))
                        {
                            auto& field_component = field_components[j];
                            dofs_to_constrain_by_field_component[j].push_back(
                              field_dof_map.at_node(i)[field_component]);
                        }
                    }
                }
            }
        }
    }
    else
    {
        throw _UTIL_MSG_EXCEPT("BCApplicationType " + to_string(application_type) +
                               " is not implemented for this function.");
    }

    return dofs_to_constrain_by_field_component;
}

void PeriodicBC::initialize(const finite_element::FiniteElementMesh& fe_mesh,
                            const mesh::NodeToElementMap& node_to_element_map,
                            const dof_manager::DofMap& dof_map,
                            const dof_manager::RankDofIndex& n_components_for_field,
                            const double& scale_factor,
                            const bool& is_incremental,
                            sparse_linear_algebra::ConstraintMatrix& constraint_matrix)
{

    const auto& rank_to_global_map = dof_map.rank_to_global_map;

    std::vector<unsigned short> field_components(n_components_for_field);
    std::iota(field_components.begin(), field_components.end(), 0);

    auto primary_dofs_by_field_component =
      find_dofs_to_constrain_by_field_component(primary_node_group_name_,
                                                application_type_,
                                                field_components,
                                                fe_mesh,
                                                node_to_element_map,
                                                dof_map);

    auto secondary_dofs_by_field_component =
      find_dofs_to_constrain_by_field_component(secondary_node_group_name_,
                                                application_type_,
                                                field_components,
                                                fe_mesh,
                                                node_to_element_map,
                                                dof_map);

    dof_manager::RankDofIndex n_primary_dofs_constrained = 0;
    for (auto& dofs : primary_dofs_by_field_component)
    {
        n_primary_dofs_constrained += dofs.size();
    }
    dof_manager::RankDofIndex n_secondary_dofs_constrained = 0;
    for (auto& dofs : secondary_dofs_by_field_component)
    {
        n_secondary_dofs_constrained += dofs.size();
    }

    _UTIL_ASSERT(n_primary_dofs_constrained != 0,
                 "Periodic boundary condition applied to " + primary_node_group_name_ + " and " +
                   secondary_node_group_name_ +
                   " did not result in any constraints. Check the boundary conditions.");

    _UTIL_ASSERT(n_primary_dofs_constrained = n_secondary_dofs_constrained,
                 "Periodic boundary condition applied to " + primary_node_group_name_ + " and " +
                   secondary_node_group_name_ +
                   " resulted in a mismatch in DoFs for the constraints. If there is an XFEM crack "
                   "on one side of the periodic boundary condition, this isn't supported yet.");

    for (auto i : utility::indices(primary_dofs_by_field_component))
    {
        for (auto j : utility::indices(primary_dofs_by_field_component[i]))
        {

            const auto& primary_global_dof =
              rank_to_global_map[primary_dofs_by_field_component[i][j]];
            const auto& secondary_global_dof =
              rank_to_global_map[secondary_dofs_by_field_component[i][j]];
            const auto& field_component = field_components[j];

            //dof_map.

            throw _UTIL_MSG_EXCEPT_C("TODO finish implementing the periodic DoFs");
            constraint_matrix.add_multi_point_constraint(primary_global_dof,
                                                         { secondary_global_dof },
                                                         { 1.0 },
                                                         0.0);
        }
    }
}
}