#ifndef __Vector3_h__
#define __Vector3_h__

#include "Common.h++"
#include "FloatUtil.h++"

//==============================================================================
// Template Vector3
//==============================================================================
template<typename Real>
struct Vector3 {
	Vector3();
	Vector3(const Vector3&);
	Vector3(const Real);
	Vector3(const Real, const Real, const Real);
	
	Real& operator[](const uint32_t);
	const Real& operator[](const uint32_t) const;

	void operator+=(const Vector3&);
	void operator-=(const Vector3&);
	void operator*=(const Real);
	void operator/=(const Real);

	Vector3 operator-() const;

	Vector3 operator+(const Vector3&) const;
	Vector3 operator-(const Vector3&) const;
	Vector3 operator*(const Real) const;
	Vector3 operator/(const Real) const;

	friend Vector3 operator*(const Real, const Vector3&);

	bool operator==(const Vector3&) const;
	bool operator!=(const Vector3&) const;

	Real* Ptr();
	const Real* Ptr() const;

	Real x, y, z;
};

//==============================================================================
// Type Definitions
//==============================================================================
typedef Vector3<char>			Vector3c;
typedef Vector3<unsigned char>	Vector3uc;
typedef Vector3<short>			Vector3s;
typedef Vector3<unsigned short>	Vector3us;
typedef Vector3<int>			Vector3i;
typedef Vector3<unsigned int>	Vector3ui;
typedef Vector3<float>			Vector3f;
typedef Vector3<double>			Vector3d;

//==============================================================================
// Global Methods
//==============================================================================
template<typename Real>
Vector3<Real> operator*(const Real, const Vector3<Real>&);

template<typename Real>
void Negate(Vector3<Real>&, const Vector3<Real>&);

template<typename Real>
void Add(Vector3<Real>&, const Vector3<Real>&, const Vector3<Real>&);

template<typename Real>
void Subtract(Vector3<Real>&, const Vector3<Real>&, const Vector3<Real>&);

template<typename Real>
void Scale(Vector3<Real>&, const Vector3<Real>&, const Real);

template<typename Real>
void ScaleAdd(Vector3<Real>&,
		   	  const Vector3<Real>&,
			  const Vector3<Real>&,
			  const Real);

template<typename Real>
bool Compare(const Vector3<Real>&, const Vector3<Real>&);

template<typename Real>
Real Length(const Vector3<Real>&);

template<typename Real>
Real LengthSquared(const Vector3<Real>&);

template<typename Real>
Real Dot(const Vector3<Real>&, const Vector3<Real>&);

template<typename Real>
void Cross(Vector3<Real>&, const Vector3<Real>&, const Vector3<Real>&);

template<typename Real>
void Normalize(Vector3<Real>&);

template<typename Real>
void Lerp(Vector3<Real>&,
		  const Vector3<Real>&,
		  const Vector3<Real>&,
		  const Real);

template<typename Real>
void Hermite(Vector3<Real>&,
			 const Vector3<Real>&,
			 const Vector3<Real>&,
			 const Vector3<Real>&,
			 const Vector3<Real>&,
			 const Real);

template<typename Real>
void CatmullRom(Vector3<Real>&,
				const Vector3<Real>&,
				const Vector3<Real>&,
				const Vector3<Real>&,
				const Vector3<Real>&,
				const Real);

template<typename Real>
void Barycentric(Vector3<Real>&,
				 const Vector3<Real>&,
				 const Vector3<Real>&,
				 const Vector3<Real>&,
				 const Real,
				 const Real);

//==============================================================================
// Template Implementation
//==============================================================================
#include "Vector3.inl"

#endif
