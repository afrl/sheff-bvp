#pragma once

#include "fea_analysis/LinearAnalysis.h"

#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/Vector.h"
#include "sparse_linear_algebra/SparseSolver.h"

#include "fea_crack/EnrichmentSurface.h"
#include "fea_crack/EnrichmentSurfaceManager.h"

#include "finite_element/Factory.h"

#include <memory>
#include <vector>

namespace fea_analysis
{
    /// For this type of surface insertion, the model for each basal fields will be unchanged,
    /// but the surface incrementally evolves based on the growth and direction criteria.
    struct SurfaceInitiationAtPoint
    {
        simple_math::Point initiation_point;
        constitutive::EnrichmentProperties& properties;
        constitutive::GrowthCriterion& growth_criterion;
        constitutive::DirectionCriterion& direction_criterion;
    };

    struct SurfaceInsertionByLevelSet
    {
        std::function<double(const simple_math::Point&)> level_set_function;
        constitutive::EnrichmentProperties& properties;
        constitutive::FieldProperties* properties_for_fields_on_neg_side;
    };

    /// A linear elastic analysis
    class FractureAnalysis : public LinearElasticAnalysis
    {
    public:
        using LevelSetFunction = std::function<double(const simple_math::Point&)>;

        FractureAnalysis();
        virtual ~FractureAnalysis() {}

        /// Creates all of the objects needed before solve_time_steps can be called.
        /// Also sets the dof map phase to enrichment so any insert_* methods can be called.
        virtual void setup_analysis() override;

        /// Insert precracks by initiating a new crack at the given point and growing the surface
        /// via the direction criterion and growth criterion that were assigned to the elements.
        void insert_precracks(
            const std::vector<SurfaceInitiationAtPoint>& initations);
        /// Insert precracks along implicit functions that define the surface (where the functions
        /// equals 0).  This method will insert surfaces wherever the function is 0, regardless of
        /// other surfaces.  For example, if the first function defines a sphere and the second
        /// function defines a plane through the sphere, a crack will be inserted along the plane
        /// both inside and outside the sphere.
        void insert_precracks_by_function(
            std::vector<SurfaceInsertionByLevelSet>& insertions);
        /// Insert precracks along implicit functions that define the surface (where the functions
        /// equals 0).  This method differs from insert_precracks_by_function in that it only
        /// inserts a surface where the sign of previously inserted surfaces (in the vector) are
        /// positive.  For example, if the first function defines a sphere with the normal pointing
        /// outward and the second function defines a plane intersecting the sphere, then a surface
        /// will only be inserted outside the sphere.
        void insert_precracks_by_function_outside_previous(
            const std::vector<SurfaceInsertionByLevelSet>& insertions);
        /// TODO document
        void insert_precracks_by_function_with_boolean(
            const std::vector<SurfaceInsertionByLevelSet>& insertions);

		/// Loops over all enrichment surfaces and tries to grow them again.  Presumably, the propagation
		/// criteria for the surface could have changed after the insert function was called, and this
		/// method will try to grow them again.
		void try_to_grow_precracks();

    protected:
        std::unique_ptr<fea_crack::EnrichmentSurfaceManager> enrichment_surfaces_;

        /// This method should include any tasks that need to be done at the beginning of
        /// FeaAnalysis::solve_time_steps().  For this class, it will change the dof map phase 
        /// from enrichment to a static phase suitable for solving (renumbering any new DoF's).
        virtual void pre_solve_tasks() override;

        virtual void solve_time_step() override;
        virtual void write_output() override;

        void solve_equillibrium_step();
        bool try_to_initiate_cracks();
        bool try_to_grow_cracks();
    };
}
