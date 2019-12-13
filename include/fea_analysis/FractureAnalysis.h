#pragma once

#include "fea_analysis/LinearAnalysis.h"

#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/Vector.h"
#include "sparse_linear_algebra/SparseSolver.h"

#include "fea_crack/EnrichmentSurface.h"
#include "fea_crack/EnrichmentSurfaceManager.h"

#include <memory>
#include <vector>

namespace fea_analysis
{
    /// A linear elastic analysis
    class FractureAnalysis : public LinearElasticAnalysis
    {
    public:
        using LevelSetFunction = std::function<const double(const simple_math::Point&)>;

        FractureAnalysis();
        virtual ~FractureAnalysis() {}

        /// Creates all of the objects needed before solve_time_steps can be called
        virtual void setup_analysis() override;

        /// Insert precracks by initiating a new crack at the given point and growing the surface
        /// via the direction criterion and growth criterion that were assigned to the elements.
        void insert_precracks(const std::vector<simple_math::Point>& initiation_points);
        /// Insert precracks along implicit functions that define the surface (where the functions
        /// equals 0).  This method will insert surfaces wherever the function is 0, regardless of
        /// other surfaces.  For example, if the first function defines a sphere and the second
        /// function defines a plane through the sphere, a crack will be inserted along the plane
        /// both inside and outside the sphere.
        void insert_precracks_by_function(
            const std::vector<LevelSetFunction>& level_set_functions);
        /// Insert precracks along implicit functions that define the surface (where the functions
        /// equals 0).  This method differs from insert_precracks_by_function in that it only
        /// inserts a surface where the sign of previously inserted surfaces (in the vector) are
        /// positive.  For example, if the first function defines a sphere with the normal pointing
        /// outward and the second function defines a plane intersecting the sphere, then a surface
        /// will only be inserted outside the sphere.
        void insert_precracks_by_function_outside_previous(
            const std::vector<LevelSetFunction>& level_set_functions);
        /// TODO document
        void insert_precracks_by_function_with_boolean(
            const std::vector<LevelSetFunction>& level_set_functions);

    protected:
        fea_crack::ElementAdjacencyList adjacency_list_;
        std::unique_ptr<fea_crack::EnrichmentSurfaceManager> enrichment_surfaces_;

        virtual void solve_time_step() override;
        virtual void write_output() override;

        void solve_equillibrium_step();
        bool try_to_initiate_cracks();
        bool try_to_grow_cracks();
    };
}