#ifndef __Coordinates_h__
#define __Coordinates_h__

#include "FloatUtil.hpp"
#include "Vector2.hpp"
#include "Vector3.hpp"

//==============================================================================
// Template Cartesian2
//==============================================================================
template<typename Real>
struct Cartesian2 : public Vector2<Real>
{
    Cartesian2()
        : Vector2<Real>()
    {
    }

    Cartesian2(const Real x0, const Real y0)
        : Vector2<Real>(x0, y0)
    {
    }

    Cartesian2(const Cartesian2<Real>& c)
        : Vector2<Real>(c.x, c.y)
    {
    }
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Cartesian2<float> Cartesian2f;
typedef Cartesian2<double> Cartesian2d;

//==============================================================================
// Template Cartesian3
//==============================================================================
template<typename Real>
struct Cartesian3 : public Vector3<Real>
{
    Cartesian3()
        : Vector3<Real>()
    {
    }

    Cartesian3(const Real x0, const Real y0, const Real z0)
        : Vector3<Real>(x0, y0, z0)
    {
    }

    Cartesian3(const Cartesian3<Real>& c)
        : Vector3<Real>(c.x, c.y, c.z)
    {
    }
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Cartesian3<float> Cartesian3f;
typedef Cartesian3<double> Cartesian3d;

//==============================================================================
// Template Polar
//==============================================================================
template<typename Real>
struct Polar
{
    Polar()
        : r(Real(0)), theta(Real(0))
    {
    }

    Polar(const float r0, const float theta0)
        : r(r0), theta(theta0)
    {
    }

    Polar(const Polar<Real>& p)
        : r(p.r), theta(p.theta)
    {
    }

    Real r, theta;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Polar<float> Polarf;
typedef Polar<double> Polard;

//==============================================================================
// Template Cylindrical
//==============================================================================
template<typename Real>
struct Cylindrical
{
    Cylindrical()
        : r(Real(0)), theta(Real(0)), z(Real(0))
    {
    }

    Cylindrical(const float r0, const float theta0, const float z0)
        : r(r0), theta(theta0), z(z0)
    {
    }

    Cylindrical(const Cylindrical<Real>& c)
        :  r(c.r), theta(c.theta), z(c.z)
    {
    }

    Real r, theta, z;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Cylindrical<float> Cylindricalf;
typedef Cylindrical<double> Cylindricald;

//==============================================================================
// Template Spherical
//==============================================================================
template<typename Real>
struct Spherical
{
    Spherical()
        : rho(Real(0)), theta(Real(0)), phi(Real(0))
    {
    }

    Spherical(const float rho0, const float theta0, const float phi0)
        : rho(rho0), theta(theta0), phi(phi0)
    {
    }

    Spherical(const Spherical<Real>& s)
        : rho(s.rho), theta(s.theta), phi(s.phi)
    {
    }

    Real rho, theta, phi;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Spherical<float> Sphericalf;
typedef Spherical<double> Sphericald;

//==============================================================================
// Conversion Methods
//==============================================================================
template<typename Real>
void PolarToCartesian(Cartesian2<Real>& c, const Polar<Real>& p);

template<typename Real>
void CartesianToPolar(Polar<Real>& p, const Cartesian2<Real>& c);

template<typename Real>
void CartesianToCylindrical(Cylindrical<Real>& cylindrical,
    const Cartesian3<Real>& cartesian);

template<typename Real>
void CylindricalToCartesian(Cartesian3<Real>& cartesian,
    const Cylindrical<Real>& cylindrical);

//==============================================================================
// Inline Implementation
//==============================================================================
#include "Coordinates.inl"

#endif
