#pragma once

#include "Defines.h"
#include "dof_manager/DofMap.h"

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
            const dof_manager::LocalIndex& n_components_for_field,
			sparse_linear_algebra::ConstraintMatrix& constraint_matrix) = 0;
	};

	struct DirchletBCOpts
	{
		DirchletBCOpts() : field_component(0), value(0) {}
		DirchletBCOpts(unsigned short field_component_a, double value_a)
			: field_component(field_component_a), value(value_a) {}

		/// Component of the field variable, default 0
		unsigned short field_component;
		/// Value to assign to the component of the field variable, default 0
		double value;
	};

	class DirchletBC : public BoundaryCondition
	{
	public:
		DirchletBC(NodeSelection selection, DirchletBCOpts opts)
			: selection_(selection), opts_(opts)
		{
		}
		virtual ~DirchletBC() {}

	protected:
		NodeSelection selection_;
		DirchletBCOpts opts_;


		/// Called by FeaAnalysis::add_boundary_condition
		virtual void update_constraints() override
		{

		}

		virtual void initialize(const Mesh& mesh,
			const dof_manager::DofMap& dof_map,
            const dof_manager::LocalIndex& n_components_for_field,
			sparse_linear_algebra::ConstraintMatrix& constraint_matrix) override
		{
			selection_.run_selection(mesh);

			auto node_dof_map = dof_map.get_array_of_node_maps();
            auto n_node_dofs = dof_map.get_array_of_node_n_dofs();
			for (auto& n_ptr : selection_.get_selection())
			{
                const auto n_dofs_for_node = n_node_dofs[n_ptr->get_local_number()];
                const auto n_dofs_to_set = n_dofs_for_node / n_components_for_field;

                for (auto i : utility::range_from_zero(n_dofs_to_set))
                {
                    auto offset = i * n_components_for_field + opts_.field_component;
                    auto global_dof = node_dof_map[n_ptr->get_local_number()][offset];

                    auto c = std::make_unique<sparse_linear_algebra::DirchletConstraint>(global_dof);
                    c->set_value(opts_.value);
                    constraint_matrix.add_constraint(std::move(c));
                }
				
			}
		}
	};

	// TODO add FunctionalDirchletBC and FunctionalDirchletBCOpts
}