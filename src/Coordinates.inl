//------------------------------------------------------------------------------
// PolarToCartesian
//------------------------------------------------------------------------------
template<typename Real>
inline void
PolarToCartesian(Cartesian2<Real>& c, const Polar<Real>& p)
{
	c.x = p.r * FloatUtil<Real>::Cos(p.theta);
	c.y = p.r * FloatUtil<Real>::Sin(p.theta);
}

//------------------------------------------------------------------------------
// CartesianToPolar
//------------------------------------------------------------------------------
template<typename Real>
inline void
CartesianToPolar(Polar<Real>& p, const Cartesian2<Real>& c)
{
	p.r = FloatUtil<Real>::Sqrt(c.x * c.x + c.y * c.y);
	p.theta = FloatUtil<Real>::ArcTan2(c.x, c.y);
}

//------------------------------------------------------------------------------
// CartesianToCylindrical
//------------------------------------------------------------------------------
template<typename Real>
inline void
CartesianToCylindrical(Cylindrical<Real>& cylindrical,
					   const Cartesian3<Real>& cartesian)
{
	cylindrical.r = FloatUtil<Real>::Sqrt(cartesian.x * cartesian.x +
										  cartesian.y * cartesian.y);
	cylindrical.theta = FloatUtil<Real>::ArcTan2(cartesian.x, cartesian.y);
	cylindrical.z = cartesian.z;
}

//------------------------------------------------------------------------------
// CylindricalToCartesian
//------------------------------------------------------------------------------
template<typename Real>
inline void
CylindricalToCartesian(Cartesian3<Real>& cartesian,
					   const Cylindrical<Real>& cylindrical)
{
	cartesian.x = cylindrical.r * FloatUtil<Real>::Cos(cylindrical.theta);
	cartesian.y = cylindrical.r * FloatUtil<Real>::Sin(cylindrical.theta);
	cartesian.z = cylindrical.z;
}