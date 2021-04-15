#include "../UnitTestHelpers.h"
#include "bvp/bvp/FractureBVP.h"
#include "constitutive/enrichment/deformation/bonded/BondedInterfaceProperties.h"
#include "constitutive/enrichment/deformation/bonded/InfStrainBondedEnrichmentModel.h"
#include "constitutive/enrichment/deformation/traction_free/InfStrainTractionFreeEnrichmentModel.h"
#include "constitutive/enrichment/deformation/traction_free/TractionFreeInterfaceProperties.h"
#include "constitutive/enrichment/direction_criteria/MaxPrincipalStressDirectionCriterion.h"
#include "constitutive/enrichment/direction_criteria/PositionFunctionDirectionCriterion.h"
#include "constitutive/enrichment/growth_criteria/AlwaysGrowGrowthCriterion.h"
#include "constitutive/field/deformation/linear_elasticity/LinearElasticProperties.h"
#include "finite_element/FiniteElement.h"
#include "finite_element/FiniteElementMesh.h"
#include "finite_element/ModelSelectionHelper.h"
#include "finite_element/element_layouts/HierarchicalFieldLayout.h"
#include "finite_element/element_layouts/SingleFieldLayout.h"
#include "finite_element/enrichment_models/deformation/DeformationEnrichmentModel.h"
#include "finite_element/enrichment_models/deformation/TractionFreeEnrichmentModel.h"
#include "finite_element/hhe_primitives/PrimitiveAlgorithms.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/UnitElements.h"
#include "simple_math/quadrature/Quadrature.h"
#include "utility/Profile.h"

using namespace bvp;
using namespace finite_element;

double compute_point_to_plane_distance(simple_math::Point origin,
                                       simple_math::Point normal,
                                       simple_math::Point point)
{
    auto o = simple_math::point_as_vector(origin);
    auto n = simple_math::point_as_vector(normal);
    auto x = simple_math::point_as_vector(point);
    return n.dot(x - o) / n.norm();
}

simple_math::GeneralVector compute_nodal_signed_distances(const mesh::ElementAccessor& element,
                                                          simple_math::Point origin,
                                                          simple_math::Vector3D normal)
{
    simple_math::GeneralVector nodal_sign_dist(element.connectivity().size());
    for (auto i : utility::range(element.connectivity().size()))
    {
        nodal_sign_dist(i) = compute_point_to_plane_distance(origin,
                                                             simple_math::to_point(normal),
                                                             element.node(i).point());
    }
    return nodal_sign_dist;
}

double sum_quad_weights(const ElementLayout& layout, const LeafIndex basal_field_index)
{
    auto& quadrature_data_0 =
      *layout.get_basal_field(basal_field_index).model_data.adaptive_quadrature;
    double weight_sum = 0.0;
    for (auto i : utility::range(quadrature_data_0.n_quad_points()))
    {
        weight_sum += quadrature_data_0.get_quad_weight(i);
    }
    return weight_sum;
}

double sum_enrichment_quad_weights(const ElementLayout& layout, const LeafIndex enrichment_index)
{
    double weight_sum = 0.0;
    for (auto& patch_data : layout.get_enrichment(enrichment_index).model_data.patch_data)
    {
        for (auto i : utility::range(patch_data.adaptive_quadrature.n_quad_points()))
        {
            weight_sum += patch_data.adaptive_quadrature.get_quad_weight(i);
        }
    }
    return weight_sum;
}


struct ParametricEntry
{
    simple_math::Vector3D position;
    simple_math::Vector3D normal;
};

int main(int argc, char** argv)
{
    fea_io::MPIUtils::init(&argc, &argv);
    utility::MPIUtils::init(&argc, &argv);
    sparse_linear_algebra::LinAlgFactory::init(
      sparse_linear_algebra::LinAlgFactory::LinAlgBackend::PETSC,
      &argc,
      &argv);

    bvp::BVPGlobalOptions::get().save_global_matrices = true;

    utility::Profile::set_quiet(true);

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Input parameters

    // Parameters to control what enrichment surface orientations and positijons are sampled in the
    // element:
    int N = 7;                      // # of points along a uniform grid edge
    int M = 7;                      // controls # of orientations sampled (\sum_{i=1}^{M + 1}{i})
    assert(M > 1);                  // The sampling algorithm results in NaN values if M <= 1
    int n_singular_eigenvalues = 6; // For the bonded problem, there is 6, but for the open crack
                                    // problem there is 12 singular values.
    // Elasticity parameters
    double E  = 1;
    double nu = 0.3;

    // Solver/PC combinations to try
    std::vector<std::pair<sparse_linear_algebra::PreconditionerCategories,
                          sparse_linear_algebra::SolverCategories>>
      solvers_to_test = { { sparse_linear_algebra::PreconditionerCategories::JACOBI,
                            sparse_linear_algebra::SolverCategories::CG },
                          { sparse_linear_algebra::PreconditionerCategories::ILU,
                            sparse_linear_algebra::SolverCategories::CG },
                          { sparse_linear_algebra::PreconditionerCategories::JACOBI,
                            sparse_linear_algebra::SolverCategories::GMRES },
                          { sparse_linear_algebra::PreconditionerCategories::ILU,
                            sparse_linear_algebra::SolverCategories::GMRES } };


    //////////////////////////////////////////////////////////////////////////////////////////////
    // Setup properties, criteria, and crack samples

    auto& field_props = constitutive::FieldPropertiesCollection::singleton();
    constitutive::LinearElasticProperties props;
    props.set_isotropic_E_nu(E, nu);
    auto props_ref = field_props.add("isotropic", std::move(props));

    auto& enrichment_prop_collection = constitutive::EnrichmentPropertiesCollection::singleton();
    auto enrichment_bond_props_ref =
      enrichment_prop_collection.add("bonded",
                                     constitutive::BondedInterfaceProperties::create(1e10));

    auto& growth_criteria_collection = constitutive::GrowthCriterionCollection::singleton();
    auto growth_criterion_ref =
      growth_criteria_collection.add("always_grow", constitutive::AlwaysGrowGrowthCriterion());

    auto& direction_criteria_collection = constitutive::DirectionCriterionCollection::singleton();
    auto direction_criterion_ref =
      direction_criteria_collection.add("max_prin_stress",
                                        constitutive::MaxPrincipalStressDirectionCriterion());

    auto& model_selection = finite_element::ModelSelectionHelper::get();

    std::unique_ptr<finite_element::ElementWorkspace> workspace_ptr =
      std::make_unique<finite_element::ElementWorkspace>();
    simple_math::GeneralVector zero_solution_vector;


    auto position_coords = simple_math::linspace(-0.99, 0.99, N);

    auto C = sin(EIGEN_PI / 4.0);

    std::vector<ParametricEntry> samples;
    for (auto i : utility::range(N + 1))
    {
        for (auto j : utility::range(N + 1))
        {
            for (auto k : utility::range(N + 1))
            {
                // std::cout << samples.size() << " (" << position_coords[i] << ", "
                //          << position_coords[j] << ", " << position_coords[k] << ")" << std::endl;
                //
                // Create a relatively evenly spaced distribution of orientations over an 1/8th of a
                // sphere Note: I sample theta and phi and then transform the points to improve the
                // evenness
                for (auto latitude : utility::range(M))
                {
                    auto phi = EIGEN_PI / 4.0 * (1.0 - (double)latitude / (double)(M - 1));

                    for (auto longitude : utility::range(M - latitude))
                    {
                        auto theta =
                          M - latitude == 1
                            ? 0.0
                            : EIGEN_PI / 2.0 * ((double)longitude / (double)(M - latitude - 1));

                        auto u      = sin(phi) * cos(theta) / C;
                        auto v      = sin(phi) * sin(theta) / C;
                        auto w      = (cos(phi) - C) / (1.0 - C);
                        auto normal = simple_math::Vector3D(u, v, w);
                        // std::cout << normal(0) << " " << normal(1) << " " << normal(2) <<
                        // std::endl;
                        samples.push_back({ simple_math::Vector3D(position_coords[i],
                                                                  position_coords[j],
                                                                  position_coords[k]),
                                            normal });
                    }
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Solve the eigenvalue problem for each sample (without boundary conditions)

    std::cout << "Finding the eigenvalues of each sample." << std::endl;
    std::size_t output_mark = samples.size() / 10;
    std::vector<std::pair<int, double>> samples_min_eigenvalues;
    for (auto i : utility::indices(samples))
    {
        auto& sample = samples[i];
        if (i % output_mark == 0)
        {
            std::cout << "   " << ((double)i / (double)samples.size()) * 100 << "% complete (at "
                      << i << ")\n";
        }

        auto& normal = sample.normal;

        FiniteElementMesh fe_mesh(mesh::get_unit_linear_hex());
        for (auto& finite_element : fe_mesh.finite_elements())
        {
            model_selection.assign_field_models(finite_element, "InfinitesimalStrain", props_ref);
            FiniteElement::allocate_isvs(finite_element.layout);
        }

        auto surface_model =
          model_selection.create_enrichment_surface_model("InfinitesimalStrain",
                                                          enrichment_bond_props_ref,
                                                          growth_criterion_ref,
                                                          direction_criterion_ref);


        auto& element_0 = fe_mesh.mesh.elements()[0];
        simple_math::GeneralVector nodal_sign_dist =
          compute_nodal_signed_distances(element_0, simple_math::to_point(sample.position), normal);

        try
        {

            Enricher::enrich_element(element_0,
                                     fe_mesh.layouts[0],
                                     0, // leaf to enrich
                                     0, // surface id
                                     nodal_sign_dist,
                                     std::vector<EnrichmentCompatibility>(),
                                     surface_model.enrichment_model,
                                     surface_model.pointwise_enrichment_model_data);


            zero_solution_vector.resize(3 * element_0.connectivity().size() *
                                        fe_mesh.layouts[0]->get_n_basal_fields());
            zero_solution_vector.setZero();

            auto params = FiniteElement::make_params(element_0,
                                                     zero_solution_vector,
                                                     *fe_mesh.layouts[0],
                                                     *workspace_ptr);
            workspace_ptr->lock();
            auto element_mat_vecs_0 =
              FiniteElement::calculate(params,
                                       utility::mask(FEMatricesAndVectors::Flags::LHSMatrix) |
                                         utility::mask(FEMatricesAndVectors::Flags::RHSVector));
            workspace_ptr->release();

            auto basal_field_0_volume = sum_quad_weights(*fe_mesh.layouts[0], 0);
            auto basal_field_1_volume = sum_quad_weights(*fe_mesh.layouts[0], 1);
            auto min_volume           = std::min(basal_field_0_volume, basal_field_1_volume);

            auto surface_volume = sum_enrichment_quad_weights(*fe_mesh.layouts[0], 0);

            auto& block_K   = element_mat_vecs_0.lhs_matrix;
            auto block_size = element_mat_vecs_0.lhs_matrix.get_block(0, 0).rows();
            simple_math::GeneralMatrix combined_K(block_K.get_n_rows() * block_size,
                                                  block_K.get_n_rows() * block_size);
            combined_K.setZero();
            if (block_K.get_n_rows() == 1)
            {
                combined_K = element_mat_vecs_0.lhs_matrix.get_block(0, 0);
            }
            else
            {
                // Ignore cross terms since they are 0 for an open crack
                combined_K.topLeftCorner(block_size, block_size) =
                  element_mat_vecs_0.lhs_matrix.get_block(0, 0);
                if (!element_mat_vecs_0.lhs_matrix.is_block_null(0, 1))
                {
                    combined_K.topRightCorner(block_size, block_size) =
                      element_mat_vecs_0.lhs_matrix.get_block(0, 1);
                }
                if (!element_mat_vecs_0.lhs_matrix.is_block_null(1, 0))
                {
                    combined_K.topRightCorner(block_size, block_size) =
                      element_mat_vecs_0.lhs_matrix.get_block(1, 0);
                }
                combined_K.bottomRightCorner(block_size, block_size) =
                  element_mat_vecs_0.lhs_matrix.get_block(1, 1);
            }

            Eigen::EigenSolver<simple_math::GeneralMatrix> eigensolver;
            eigensolver.compute(combined_K);
            auto eigenvalues = eigensolver.eigenvalues();

            std::vector<double> sorted_eigenvalues;
            for (auto i : utility::range(eigenvalues.size()))
            {
                sorted_eigenvalues.push_back(eigenvalues(i).real());
            }
            std::sort(sorted_eigenvalues.begin(), sorted_eigenvalues.end());

            samples_min_eigenvalues.push_back({ i, sorted_eigenvalues[n_singular_eigenvalues] });
        }
        catch (...)
        {
            std::cout << "Sample " << i << " failed: " << std::endl
                      << "  position = (" << sample.position(0) << ", " << sample.position(1)
                      << ", " << sample.position(2) << ")" << std::endl
                      << "  normal = (" << sample.normal(0) << ", " << sample.normal(1) << ", "
                      << sample.normal(2) << ")" << std::endl;
        }
    }

    // Sort the sample by the smallest non-singular eigenvalue (rigid body modes are expected since
    // no boundary conditions are applied)
    std::sort(samples_min_eigenvalues.begin(), samples_min_eigenvalues.end());

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Solve the boundary value problem for each sample (sorted by the smallest non-singular
    // eigenvalue) noting whether or not the solver succeeds.

    simple_math::Vector3D curr_normal;
    auto updatable_direction_criterion_ref = direction_criteria_collection.add(
      "updatable_direction_criterion",
      constitutive::PositionFunctionDirectionCriterion(
        [&curr_normal](const int& crack_index, const simple_math::Point& p)
          -> simple_math::Vector3D { return curr_normal; }));

    std::cout << "Testing solvers on a BVP for each sample" << std::endl;
    std::size_t count = 0;
    output_mark       = samples_min_eigenvalues.size() / 10;
    std::ofstream out("solver_results.csv");
    out << "i,x,y,z,n_x,n_y,n_z,min_nonsingular_eigenvalue,succeeded,solver_type,pc_type,n_"
           "iterations,residual_norm,time (# of microseconds)"
        << std::endl;
    for (auto& sample_min_eigenvalue : samples_min_eigenvalues)
    {
        auto& sample         = samples[sample_min_eigenvalue.first];
        auto& min_eigenvalue = sample_min_eigenvalue.second;
        curr_normal          = sample.normal;

        if (count % output_mark == 0)
        {
            std::cout << "   " << ((double)count / (double)samples_min_eigenvalues.size()) * 100.0
                      << "% complete (at " << count << ")\n";
        }

        for (auto& solver : solvers_to_test)
        {

            auto bvp = FractureBVP::get();
            BVPData bvp_data;
            auto& solver_opts                    = bvp_data.core_data.solver_opts;
            solver_opts.verbose_level            = sparse_linear_algebra::VerboseLevel::None;
            solver_opts.preffered_preconditioner = solver.first;
            solver_opts.preferred_solver         = solver.second;
            solver_opts.max_iterations           = 200;
            sparse_linear_algebra::SolverStatus solver_status;
            solver_opts.status_callback =
              [&solver_status](const sparse_linear_algebra::SolverStatus& status) {
                  solver_status = status;
              };

            bvp.set_input(std::move(mesh::get_unit_linear_hex()),
                          get_test_name(__FILE__),
                          (get_output_directory(__FILE__) / "output").string(),
                          bvp_data);
            bvp.set_field_models(finite_element::BasalFieldSelection(),
                                 "InfinitesimalStrain",
                                 "isotropic",
                                 bvp_data);
            add_unit_tension_x_bcs(bvp_data);

            auto& fe_mesh = bvp_data.fe_mesh();


            try
            {
                bvp.insert_precracks(bvp_data,
                                     { { // surface 0
                                         simple_math::to_point(sample.position),
                                         enrichment_bond_props_ref,
                                         growth_criterion_ref,
                                         updatable_direction_criterion_ref } });

                auto start = std::chrono::high_resolution_clock::now();

                bvp.solve_time_steps(bvp_data);

                auto stop = std::chrono::high_resolution_clock::now();

                out << sample_min_eigenvalue.first << "," << sample.position.x() << ","
                    << sample.position.y() << "," << sample.position.z() << "," << sample.normal.x()
                    << "," << sample.normal.y() << "," << sample.normal.z() << "," << min_eigenvalue
                    << "," << 1 << "," << solver_status.string_stats["ksp_type"] << ","
                    << solver_status.string_stats["pc_type"] << ","
                    << solver_status.integral_stats["n_iterations"] << ","
                    << solver_status.float_stats["residual_norm"] << ","
                    << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
                    << std::endl;
            }
            catch (...)
            {
                out << sample_min_eigenvalue.first << "," << sample.position.x() << ","
                    << sample.position.y() << "," << sample.position.z() << "," << sample.normal.x()
                    << "," << sample.normal.y() << "," << sample.normal.z() << "," << min_eigenvalue
                    << "," << 0 << "," << std::endl;

                std::cout << "Sample " << count << " failed: " << std::endl
                          << "  pc = " << sparse_linear_algebra::to_string(solver.first)
                          << " solver = " << sparse_linear_algebra::to_string(solver.second)
                          << std::endl
                          << "  position = (" << sample.position(0) << ", " << sample.position(1)
                          << ", " << sample.position(2) << ")" << std::endl
                          << "  normal = (" << sample.normal(0) << ", " << sample.normal(1) << ", "
                          << sample.normal(2) << ")" << std::endl;
            }
        }
        count++;
    }
    out.close();

    return 0;
}