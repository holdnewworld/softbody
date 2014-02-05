//------------------------------------------------------------------------------
// FloatUtil::Sign
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::Sign(const Real x)
{
    if (x < Real(0)) {
        return Real(-1);
    }
    if (x > Real(0)) {
        return Real(1);
    }
    return Real(0);
}

//------------------------------------------------------------------------------
// FloatUtil::Abs
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::Abs(const Real x)
{
    return x < Real(0) ? -x : x;
}

//------------------------------------------------------------------------------
// FloatUtil::Min
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::Min(const Real a, const Real b)
{
    return a < b ? a : b;
}

//------------------------------------------------------------------------------
// FloatUtil::Min
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::Min(const Real a, const Real b, const Real c)
{
    if (a < b) {
        if (a < c) {
            return a;
        } else {
            return c;
        }
    } else {
        if (b < c) {
            return b;
        } else {
            return c;
        }
    }
}

//------------------------------------------------------------------------------
// FloatUtil::Max
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::Max(const Real a, const Real b)
{
    return a < b ? b : a;
}

//------------------------------------------------------------------------------
// FloatUtil::Max
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::Max(const Real a, const Real b, const Real c)
{
    if (a > b) {
        if (a > c) {
            return a;
        } else {
            return c;
        }
    } else {
        if (b > c) {
            return b;
        } else {
            return c;
        }
    }
}


//------------------------------------------------------------------------------
// FloatUtil::NearlyZero
//------------------------------------------------------------------------------
template<typename Real>
inline bool
FloatUtil<Real>::NearlyZero(const Real x)
{
    return Abs(x) < EPSILON;
}

//------------------------------------------------------------------------------
// FloatUtil::NearlyEqual
//------------------------------------------------------------------------------
template<typename Real>
inline bool
FloatUtil<Real>::NearlyEqual(const Real x, const Real y)
{
    return Abs(x - y) < EPSILON;
}

//------------------------------------------------------------------------------
// FloatUtil::ToDegrees
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::ToDegrees(const Real radians)
{
    return RAD_TO_DEG * radians;
}

//------------------------------------------------------------------------------
// FloatUtil::ToRadians
//------------------------------------------------------------------------------
template<typename Real>
inline Real
FloatUtil<Real>::ToRadians(const Real degrees)
{
    return DEG_TO_RAD * degrees;
}

//------------------------------------------------------------------------------
// FloatUtil::Pow
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::Pow(const float x, const float y)
{
    return ::powf(x, y);
}

//------------------------------------------------------------------------------
// FloatUtil::Pow
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::Pow(const double x, const double y)
{
    return ::pow(x, y);
}

//------------------------------------------------------------------------------
// FloatUtil::Sqrt
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::Sqrt(const float x)
{
    return ::sqrtf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Sqrt
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::Sqrt(const double x)
{
    return ::sqrt(x);
}

//------------------------------------------------------------------------------
// FloatUtil::InverseSqrt
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::InverseSqrt(float x)
{
    float halfX = 0.5f * x;
    int i = *(int*)&x;
    i = 0x5f3759df - (i >> 1);
    x = *(float*)&i;
    x = x * (1.5f - halfX * x * x);
    return x;
}

//------------------------------------------------------------------------------
// FloatUtil::InverseSqrt
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::InverseSqrt(double x)
{
    return 1.0 / ::sqrt(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Lerp
//------------------------------------------------------------------------------
template<typename Real>
inline void
FloatUtil<Real>::Lerp(Real& r, const Real a, const Real b, const Real s)
{
    r = a + s * (b - a);
}

//------------------------------------------------------------------------------
// FloatUtil::Exp
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::Exp(const float x)
{
    return ::expf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Exp
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::Exp(const double x)
{
    return ::exp(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Ln
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::Ln(const float x)
{
    return ::logf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Ln
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::Ln(const double x)
{
    return ::log(x);
}

//------------------------------------------------------------------------------
// FloatUtil::Sin
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::Sin(const float radians)
{
    return ::sinf(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::Sin
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::Sin(const double radians)
{
    return ::sin(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::Cos
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::Cos(const float radians)
{
    return ::cosf(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::Cos
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::Cos(const double radians)
{
    return ::cos(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::Tan
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::Tan(const float radians)
{
    return ::tanf(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::Tan
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::Tan(const double radians)
{
    return ::tan(radians);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcSin
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::ArcSin(const float x)
{
    return ::asinf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcSin
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::ArcSin(const double x)
{
    return ::asin(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcCos
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::ArcCos(const float x)
{
    return ::acosf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcCos
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::ArcCos(const double x)
{
    return ::acos(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcTan
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::ArcTan(const float x)
{
    return ::atanf(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcTan
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::ArcTan(const double x)
{
    return ::atan(x);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcTan2
//------------------------------------------------------------------------------
template<>
inline float
FloatUtil<float>::ArcTan2(const float x, const float y)
{
    return ::atan2f(x, y);
}

//------------------------------------------------------------------------------
// FloatUtil::ArcTan2
//------------------------------------------------------------------------------
template<>
inline double
FloatUtil<double>::ArcTan2(const double x, const double y)
{
    return ::atan2(x, y);
}
