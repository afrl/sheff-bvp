#pragma once

#include "bvp/bvp/LinearElasticBVP.h"
#include "sparse_linear_algebra/SparseMatrix.h"
#include "sparse_linear_algebra/SparseSolver.h"
#include "sparse_linear_algebra/Vector.h"

#include <memory>
#include <vector>

namespace bvp
{
/// For this type of surface insertion, the model for each basal fields will be unchanged,
/// but the surface incrementally evolves based on the growth and direction criteria.
struct SurfaceInitiationAtPoint
{
    simple_math::Point initiation_point;
    constitutive::EnrichmentPropertiesCollectionRef properties;
    constitutive::GrowthCriterionCollectionRef growth_criterion;
    constitutive::DirectionCriterionCollectionRef direction_criterion;
};

struct SurfaceInsertionByLevelSet
{
    std::function<double(const simple_math::Point&)> level_set_function;
    constitutive::EnrichmentPropertiesCollectionRef properties;
    constitutive::FieldPropertiesCollectionRef properties_for_fields_on_neg_side;
};

/// A linear elastic analysis
class FractureBVP : public LinearElasticBVP
{
  public:
    using LevelSetFunction = std::function<double(const simple_math::Point&)>;

    static FractureBVP& get()
    {
        static FractureBVP singleton;
        return singleton;
    }

    virtual ~FractureBVP() {}

    /// Insert precracks by initiating a new crack at the given point and growing the surface
    /// via the direction criterion and growth criterion that were assigned to the elements.
    void insert_precracks(BVPData& bvp_data,
                          const std::vector<SurfaceInitiationAtPoint>& initations) const;

    /// Insert precracks along implicit functions that define the surface (where the functions
    /// equals 0).  This method will insert surfaces wherever the function is 0, regardless of
    /// other surfaces.  For example, if the first function defines a sphere and the second
    /// function defines a plane through the sphere, a crack will be inserted along the plane
    /// both inside and outside the sphere.
    void insert_precracks_by_function(BVPData& bvp_data,
                                      std::vector<SurfaceInsertionByLevelSet>& insertions) const;

    /// Insert precracks along implicit functions that define the surface (where the functions
    /// equals 0).  This method differs from insert_precracks_by_function in that it only
    /// inserts a surface where the sign of previously inserted surfaces (in the vector) are
    /// positive.  For example, if the first function defines a sphere with the normal pointing
    /// outward and the second function defines a plane intersecting the sphere, then a surface
    /// will only be inserted outside the sphere.
    void insert_precracks_by_function_outside_previous(
      BVPData& bvp_data,
      const std::vector<SurfaceInsertionByLevelSet>& insertions) const;

    /// TODO document
    void insert_precracks_by_function_with_boolean(
      BVPData& bvp_data,
      const std::vector<SurfaceInsertionByLevelSet>& insertions) const;

    /// Loops over all enrichment surfaces and tries to grow them again.  Presumably, the
    /// propagation criteria for the surface could have changed after the insert function was
    /// called, and this method will try to grow them again.
    bool try_to_grow_precracks(BVPData& bvp_data) const;

    virtual SolveTimeStepResult solve_next_time_step(BVPData& bvp_data) const override;

  protected:
    FractureBVP() {}

    /// @see BoundaryValueProblem::init_all_data
    virtual void init_all_data(BVPData& bvp_data) const override;

    /// @brief If enrichment has occurred, this renumbers the DoFs to satisfy the needs of the solvers
    /// and enforce continuity across rank boundaries.
    void renumber_dofs(BVPData& bvp_data) const;

    xfem::Params get_xfem_params(BVPData& bvp_data) const;
};
}
