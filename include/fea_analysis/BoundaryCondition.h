#pragma once

#include "Defines.h"
#include "dof_manager/DofMap.h"
#include "sparse_linear_algebra/ConstraintMatrix.h"
#include "finite_element/hhe_primitives/Algorithms.h"
#include "fea_crack/misc/GeneralDofUtilities.h"

namespace fea_analysis
{
	class BoundaryCondition
	{
	public:
		BoundaryCondition() {}
		~BoundaryCondition() {}

		friend class FeaAnalysis;

	protected:
		/// Called by FeaAnalysis::add_boundary_condition
		virtual void update_constraints() = 0;
		/// Called by FeaAnalysis::add_boundary_condition
		virtual void initialize(const Mesh& mesh,
			const dof_manager::DofMap& dof_map,
            const dof_manager::RankDofIndex& n_components_for_field,
			sparse_linear_algebra::ConstraintMatrix& constraint_matrix) = 0;
	};

	struct DirichletBCOpts
	{
        DirichletBCOpts() : field_component(0), value(0) {}
        DirichletBCOpts(unsigned short field_component_a, double value_a)
			: field_component(field_component_a), value(value_a) {}

		/// Component of the field variable, default 0
		unsigned short field_component;
		/// Value to assign to the component of the field variable, default 0
		double value;
	};

	class DirichletBC : public BoundaryCondition
	{
	public:
		DirichletBC(NodeSelection selection, DirichletBCOpts opts)
			: selection_(selection), opts_(opts)
		{
		}
		virtual ~DirichletBC() {}

	protected:
		NodeSelection selection_;
		DirichletBCOpts opts_;


		/// Called by FeaAnalysis::add_boundary_condition
		virtual void update_constraints() override
		{

		}

		virtual void initialize(const Mesh& mesh,
			const dof_manager::DofMap& dof_map,
            const dof_manager::RankDofIndex& n_components_for_field,
			sparse_linear_algebra::ConstraintMatrix& constraint_matrix) override
		{
			selection_.run_selection(mesh);
			auto selected_faces = flex_mesh::find_faces_from_node_set(mesh, selection_.get_selection());

			std::set<dof_manager::RankDofIndex> dofs_to_constrain;
			const auto& rank_to_global_map = dof_map.rank_to_global_map;
			for (auto& face : selected_faces)
			{
				auto& parent_element = *mesh.elements()[face.local_element_index];
				auto n_basal_fields = parent_element.get_n_basal_fields();
				fea_crack::ElementXFEMDofMap element_dof_map(dof_map, n_basal_fields, parent_element);
				const std::vector<flex_mesh::FaceInfo> face_info = flex_mesh::get_face_info(parent_element.get_type());

				for (auto basal_field_index : utility::range_from_zero(n_basal_fields))
				{
					auto field_touches_face = finite_element::does_basal_field_touch_topology(
						parent_element,
						basal_field_index,
						{ flex_mesh::TopologyType::Face, face.face_index });
					if (field_touches_face)
					{
						auto field_dof_map = element_dof_map.at_leaf(basal_field_index);
						
						// Loop over nodes in face and add rank DoF's to set
						for (auto i : face_info[face.face_index].connectivity)
						{
							dofs_to_constrain.insert(field_dof_map.at_node(i)[opts_.field_component]);
						}
					}
				}
			}

			for (auto rank_dof : dofs_to_constrain)
			{
				const auto& global_dof = rank_to_global_map[rank_dof];
				constraint_matrix.add_dirichlet_constraint(global_dof, opts_.value);
			}
		}
	};

	// TODO add FunctionalDirchletBC and FunctionalDirchletBCOpts
}