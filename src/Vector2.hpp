#ifndef __Vector2_h__
#define __Vector2_h__

#include "FloatUtil.hpp"

//==============================================================================
// Template Vector2
//==============================================================================
template<typename Real>
struct Vector2 {
    Vector2();
    Vector2(const Vector2&);
    Vector2(const Real);
    Vector2(const Real, const Real);

    Real& operator[](const int);
    const Real& operator[](const int) const;

    void operator+=(const Vector2&);
    void operator-=(const Vector2&);
    void operator*=(const Real);
    void operator/=(const Real);

    Vector2 operator-() const;

    Vector2 operator+(const Vector2&) const;
    Vector2 operator-(const Vector2&) const;
    Vector2 operator*(const Real) const;
    Vector2 operator/(const Real) const;

    friend Vector2 operator*(const Real, const Vector2&);

    bool operator==(const Vector2&) const;
    bool operator!=(const Vector2&) const;

    Real* Ptr();
    const Real* Ptr() const;

    Real x, y;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Vector2<char>           Vector2c;
typedef Vector2<unsigned char>  Vector2uc;
typedef Vector2<short>          Vector2s;
typedef Vector2<unsigned short> Vector2us;
typedef Vector2<int>            Vector2i;
typedef Vector2<unsigned int>   Vector2ui;
typedef Vector2<float>          Vector2f;
typedef Vector2<double>         Vector2d;

//==============================================================================
// Global Methods
//==============================================================================
template<typename Real>
void Negate(Vector2<Real>&, const Vector2<Real>&);

template<typename Real>
void Add(Vector2<Real>&, const Vector2<Real>&, const Vector2<Real>&);

template<typename Real>
void Subtract(Vector2<Real>&, const Vector2<Real>&, const Vector2<Real>&);

template<typename Real>
void Scale(Vector2<Real>&, const Vector2<Real>&, const Real);

template<typename Real>
void ScaleAdd(Vector2<Real>&,
              const Vector2<Real>&,
              const Vector2<Real>&,
              const Real);

template<typename Real>
bool Compare(const Vector2<Real>&, const Vector2<Real>&);

template<typename Real>
Real Length(const Vector2<Real>&);

template<typename Real>
Real InverseLength(const Vector2<Real>&);

template<typename Real>
Real LengthSquared(const Vector2<Real>&);

template<typename Real>
Real Dot(const Vector2<Real>&, const Vector2<Real>&);

template<typename Real>
void Normalize(Vector2<Real>&);

template<typename Real>
void Lerp(Vector2<Real>&,
          const Vector2<Real>&,
          const Vector2<Real>&,
          const Real);

template<typename Real>
void Hermite(Vector2<Real>&,
             const Vector2<Real>&,
             const Vector2<Real>&,
             const Vector2<Real>&,
             const Vector2<Real>&,
             const Real);

template<typename Real>
void CatmullRom(Vector2<Real>&,
                const Vector2<Real>&,
                const Vector2<Real>&,
                const Vector2<Real>&,
                const Vector2<Real>&,
                const Real);

template<typename Real>
void Barycentric(Vector2<Real>&,
                 const Vector2<Real>&,
                 const Vector2<Real>&,
                 const Vector2<Real>&,
                 const Real,
                 const Real);

//==============================================================================
// Template Implementation
//==============================================================================
#include "Vector2.inl"

#endif
