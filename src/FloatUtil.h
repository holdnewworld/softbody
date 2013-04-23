//==============================================================================
// Author: Tim Thirion
//
// This file is: FloatUtil.h
//   Created on: July 18, 2004
//  Modified on: October 18, 2004
//
//  Description:
//				Floating point constants and utility functions are defined in
//				this file. It should be noted that the class FloatUtil should
//				only use floating-point inherent C++ types as its sole template
//				argument. Some routines would work with integers, but don't
//				count on it.
//
//    Revisions:
//				1. July 17, 2004
//					FloatUtil.h created and added to the system. Constants
//					PI, DEG_TO_RAD, RAD_TO_DEG, EPSILON, and INFINITY added.
//					The methods Abs, IsZero, Compare, ToDegrees, ToRadians,
//					Sqrt, and RSqrt also added. The implementation for
//					FloatUtil<float>::RSqrt is a hack taken from the Quake III
//					source code. It is efficient although the code looks quite
//					"magical." There has been some study of it; more information
//					to follow.
//
//				2. October 18, 2004
//					Added the methods Min, Max and Pow. Changed the semantics of
//					Lerp and changed the name of IsZero to NearlyZero.
//
//==============================================================================
#ifndef __FloatUtil_h__
#define __FloatUtil_h__

#include <cmath>
#include <cassert>

const double PI = 3.141592653589793238462643383279;
const double DEG_TO_RAD = 0.017453292519943295769236907685;
const double RAD_TO_DEG = 57.29577951308232087679815481411;
const double EPSILON = 1e-3;

//==============================================================================
// class template FloatUtil
//==============================================================================
template<typename Real>
struct FloatUtil
{
	//==========================================================================
	// Public Methods
	//==========================================================================
	static Real Sign(const Real);
	static Real Abs(const Real);
	static Real Min(const Real, const Real);
	static Real Min(const Real, const Real, const Real);
	static Real Max(const Real, const Real);
	static Real Max(const Real, const Real, const Real);
	static bool NearlyZero(const Real);
	static bool NearlyEqual(const Real, const Real);
	static Real ToDegrees(const Real);
	static Real ToRadians(const Real);
	static Real Pow(const Real, const Real);
	static Real Sqrt(const Real);
	static Real InverseSqrt(Real);
	static void Lerp(Real&, const Real, const Real, const Real);

	static Real Exp(const Real);
	static Real Ln(const Real);

	static Real Sin(const Real);
	static Real Cos(const Real);
	static Real Tan(const Real);
	static Real ArcSin(const Real);
	static Real ArcCos(const Real);
	static Real ArcTan(const Real);
	static Real ArcTan2(const Real, const Real);
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef FloatUtil<float>	FloatUtilf;
typedef FloatUtil<double>	FloatUtild;

//==============================================================================
// Template Implementation
//==============================================================================
#include "FloatUtil.inl"

#endif
