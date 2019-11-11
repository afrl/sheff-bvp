#include "fea_analysis/FeaAnalysis.h"

#include "assembler/SerialAssembler.h"
#include "fea_io/meta/MetaReader.h"
#include "flex_mesh/algorithms/IOFeaIo.h"
#include "finite_element/BaseFormulation.h"

namespace fea_analysis
{
    void FeaAnalysis::set_input_meta(const std::string relative_path_to_meta)
    {
        fea_io::MetaReader reader;
        reader.open(relative_path_to_meta);
        std::string mesh_filepath = reader.get_mesh_path(0);

        mesh_ = std::make_unique<Mesh>();
        flex_mesh::FEAIOInterface<Mesh>::read_partition(*mesh_, relative_path_to_meta, 0, 0);

        output_meta_.open("output", "post-analysis");
        output_meta_.add_mesh_path(0, 0, mesh_filepath);
    }

    void FeaAnalysis::set_mesh(std::unique_ptr<Mesh> mesh)
    {
        mesh_ = std::move(mesh);

        output_meta_.open("output", "post-analysis");
        auto mesh_filepath = output_meta_.add_mesh_path(0, 0, "initial_mesh");

        flex_mesh::FEAIOInterface<Mesh>::write(*mesh_, mesh_filepath);
    }

    void FeaAnalysis::add_constitutive_relation(const std::string name, std::unique_ptr<constitutive::BaseConstitutiveRelation> relation)
    {
        constitutive_relations_[name] = std::move(relation);
    }

	void FeaAnalysis::set_elements_constitutive_relation(ElementSelection selection, const std::string constitutive_relation_name)
	{
		auto find_iter = constitutive_relations_.find(constitutive_relation_name);
		if (find_iter == constitutive_relations_.end())
		{
			throw _UTIL_MSG_EXCEPT_C("Constitutive relation \"" + constitutive_relation_name + "\" not found.");
		}

		selection.run_selection(*mesh_);

		auto constitutive_relation_ptr = find_iter->second.get();
		for (auto e_ptr : selection.get_selection())
		{
			e_ptr->set_constitutive_relation(constitutive_relation_ptr);
		}
	}

	void FeaAnalysis::add_boundary_condition(std::unique_ptr<BoundaryCondition> boundary_condition)
	{
		boundary_conditions_.push_back(std::move(boundary_condition));
	}

    void fea_analysis::FeaAnalysis::setup_analysis()
    {
        // If no time step manager has been set through the input script, it is assumed that
        // this is a static analysis.
        if (time_step_manager_ == 0)
            time_step_manager_ = std::make_unique<StaticAnalysisTimeStepManager>();

        assembler_ = std::make_unique<assembler::SerialAssembler>(*mesh_);

        set_element_formulations();

        // Setup the dof map and any global matrices/vectors.  The dof map will be
        // used to setup the containers.
        dof_map_ = std::make_unique<dof_manager::DofMap>(*mesh_);
        initialize_linear_algebra();
    }

    void FeaAnalysis::solve_time_steps()
    {
        // Time loop
        time_step_manager_->start_time_step_loop();
        do
        {
            /* TODO add restart capability */
            solve_time_step();
            write_output();
        } while (time_step_manager_->next_time_step());
    }

    void FeaAnalysis::apply_boundary_conditions()
    {
        for (auto& bc : boundary_conditions_)
        {
            // IMPORTANT NOTE: this assumes the field has # of components as dimensionality
            // This works for displacement with solid elements, but reconsider what should
            // be done to handle scalar fields, shell elements, combined fields, etc.
            // TODO figure out a robust term strategy for this
            bc->initialize(*mesh_, *dof_map_,
                utility::GlobalVariables::get_dimensionality(),
                *gl_constraint_matrix_);
        }
    }
}