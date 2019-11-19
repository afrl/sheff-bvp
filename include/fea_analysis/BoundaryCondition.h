#pragma once

#include "Defines.h"
#include "dof_manager/DofMap.h"
#include "sparse_linear_algebra/ConstraintMatrix.h"

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

			for (auto& n_ptr : selection_.get_selection())
			{
                auto node_dof_map = dof_map.get_node_map(n_ptr->get_local_number());
                const auto n_dofs_for_node = node_dof_map.size();
                const auto n_dofs_to_set = n_dofs_for_node / n_components_for_field;

                for (auto i : utility::range_from_zero(n_dofs_to_set))
                {
                    auto offset = i * n_components_for_field + opts_.field_component;
                    auto global_dof = node_dof_map[offset];
                    constraint_matrix.add_dirichlet_constraint(global_dof, opts_.value);
                }
				
			}
		}
	};

	// TODO add FunctionalDirchletBC and FunctionalDirchletBCOpts
}