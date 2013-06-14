#ifndef __Primitives_h__
#define __Primitives_h__

#include "Coordinates.h++"
#include "Vector3.h++"

//==============================================================================
// Template PointMass
//==============================================================================
template<typename Real>
struct PointMass
{
	//--------------------------------------------------------------------------
	// PointMass
	//--------------------------------------------------------------------------
	PointMass()
	: mass(0),
	  position(),
	  velocity()
	{
	}

	//--------------------------------------------------------------------------
	// PointMass
	//--------------------------------------------------------------------------
	PointMass(const Real _mass,
			  const Cartesian3<Real> _position,
			  const Vector3<Real> _normal,
			  const Vector3<Real> _velocity,
			  const Vector3<Real> _force)
	: mass(_mass),
	  position(_position),
	  normal(_normal),
	  velocity(_velocity),
	  force(_force)
	{
	}

	Real				mass;
	Cartesian3<Real>	position;
	Vector3<Real>		normal;
	Vector3<Real>		velocity;
	Vector3<Real>		force;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef PointMass<float>		PointMassf;
typedef PointMass<double>		PointMassd;

//==============================================================================
// Template SpringEdge
//==============================================================================
template<typename Real>
struct SpringEdge
{
	SpringEdge()
	: springConstant(Real(0)),
	  restLength(Real(0))
	{
		vertices[0] = vertices[1] = 0;
	}

	uint32_t	vertices[2];
	Real		springConstant;
	Real		restLength;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef SpringEdge<float>	SpringEdgef;
typedef SpringEdge<double>	SpringEdged;

//==============================================================================
// Template Face
//==============================================================================
template<typename Real>
struct Face
{
	Face()
	: normal()
	{
		vertices[0] = vertices[1] = vertices[2] = 0;
	}

	uint32_t	vertices[3];
	Vector3f	normal;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Face<float>		Facef;
typedef Face<double>	Faced;

#endif
