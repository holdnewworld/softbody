#ifndef __Plane_h__
#define __Plane_h__

#include "Coordinates.h"
#include "Vector3.h"

//==============================================================================
// Enum PlaneSide
//==============================================================================
enum PlaneSide
{
	PLANE_SIDE_NEGATIVE,
	PLANE_SIDE_ON,
	PLANE_SIDE_POSITIVE
};

//==============================================================================
// Template Plane
//==============================================================================
template<typename Real>
struct Plane
{
	Plane()
	{
	}

	Plane(const Vector3<Real>& _normal, const Cartesian3<Real> p0)
	: normal(_normal),
	  p(p0)
	{
	}

	PlaneSide PointOnPlaneSide(const Cartesian3<Real>& point)
	{
		Real test = Dot(normal, point) - Dot(normal, p);
		if (FloatUtil<Real>::NearlyZero(test)) {
			return PLANE_SIDE_ON;
		}
		if (test > Real(0)) {
			return PLANE_SIDE_POSITIVE;
		}
		return PLANE_SIDE_NEGATIVE;
	}

	Vector3<Real> normal;
	Cartesian3<Real> p;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Plane<float>	Planef;
typedef Plane<double>	Planed;

#endif