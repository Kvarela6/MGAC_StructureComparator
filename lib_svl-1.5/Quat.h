#ifndef __Quat__
#define __Quat__
//#include "svl/Vec3.h"
//#include "svl/Mat3.h"

class Quat 
{
public:
	// constructors
	Quat();
	Quat(Real q0, Real q1, Real q2, Real q3); // [q0,(q1,q2,q3)]
	Quat (const Vec3 &axis, Real angle);
	Quat (const Mat3 &m);
	
	Int	Elts() const { return (4); };
	
    Real        &operator [] (Int i);
    const Real  &operator [] (Int i) const;
	
    // Assignment operators

    Quat        &operator =  (const Quat &a);
    Quat        &operator += (const Quat &a);
    Quat        &operator -= (const Quat &a);
    Quat        &operator *= (const Quat &a);
    Quat        &operator *= (Real s);
    Quat        &operator /= (Real s);

    // Arithmetic operators

    Quat        operator + (const Quat &a) const;   // v + a
    Quat        operator - (const Quat &a) const;   // v - a
    Quat        operator - () const;                // -v
    Quat        operator * (const Quat &a) const;   // v * a (vx * ax, ...)
    Quat        operator * (Real s) const;          // v * s
    Quat        operator / (Real s) const;          // v / s


    Quat        &Normalise();                       // normalise vector

protected:
	Real elt[4];
};

inline Quat     operator * (Real s, const Quat &v); // Left mult. by s
inline Real     dot(const Quat &a, const Quat &b);  // v . a
inline Real     len(const Quat &v);                 // || v ||
inline Real     sqrlen(const Quat &v);              // v . v
inline Quat     norm(const Quat &v);                // v / || v ||
inline Void     normalise(Quat &v);                 // v = norm(v)
inline Quat		slerp(const Quat &q1, const Quat &q2, Real t);
inline Quat     conjugate(const Quat &q);

Mat3     Rot3(const Quat &q);
Mat4     HRot4(const Quat &q);

			
std::ostream &operator << (std::ostream &s, const Quat &v);
std::istream &operator >> (std::istream &s, Quat &v);

inline Real &Quat::operator [] (Int i)
{
    CheckRange(i, 0, 4, "(Quat::[i]) index out of range");
    return(elt[i]);
}

inline const Real &Quat::operator [] (Int i) const
{
    CheckRange(i, 0, 4, "(Quat::[i]) index out of range");
    return(elt[i]);
}

inline Quat::Quat()
{
}

inline Quat::Quat(Real q0, Real q1, Real q2, Real q3)
{
	elt[0] = q0;
	elt[1] = q1;
	elt[2] = q2;
	elt[3] = q3;
}

inline Quat::Quat(const Vec3 &axis, Real angle)
{
	Vec3 n = norm(axis);
	Real sinhalf = sin(angle/2);
	elt[1] = sinhalf*n[0];
	elt[2] = sinhalf*n[1];
	elt[3] = sinhalf*n[2];

	elt[0] = cos(angle/2);
}


inline Quat &Quat::operator = (const Quat &v)
{
    elt[0] = v[0];
    elt[1] = v[1];
    elt[2] = v[2];
    elt[3] = v[3];

    return(SELF);
}

inline Quat &Quat::operator += (const Quat &v)
{
    elt[0] += v[0];
    elt[1] += v[1];
    elt[2] += v[2];
    elt[3] += v[3];

    return(SELF);
}

inline Quat &Quat::operator -= (const Quat &v)
{
    elt[0] -= v[0];
    elt[1] -= v[1];
    elt[2] -= v[2];
    elt[3] -= v[3];

    return(SELF);
}

inline Quat &Quat::operator *= (const Quat &v)
{
	Quat tmp(elt[0],elt[1],elt[2],elt[3]);
	tmp = tmp * v;

    elt[0] = tmp[0];
    elt[1] = tmp[1];
    elt[2] = tmp[2];
    elt[3] = tmp[3];
	
    return(SELF);
}

inline Quat &Quat::operator *= (Real s)
{
    elt[0] *= s;
    elt[1] *= s;
    elt[2] *= s;
    elt[3] *= s;

    return(SELF);
}

inline Quat &Quat::operator /= (Real s)
{
    elt[0] /= s;
    elt[1] /= s;
    elt[2] /= s;
    elt[3] /= s;

    return(SELF);
}


inline Quat Quat::operator + (const Quat &a) const
{
    Quat result;

    result[0] = elt[0] + a[0];
    result[1] = elt[1] + a[1];
    result[2] = elt[2] + a[2];
    result[3] = elt[3] + a[3];

    return(result);
}

inline Quat Quat::operator - (const Quat &a) const
{
    Quat result;

    result[0] = elt[0] - a[0];
    result[1] = elt[1] - a[1];
    result[2] = elt[2] - a[2];
    result[3] = elt[3] - a[3];

    return(result);
}

inline Quat Quat::operator - () const
{
    Quat result;

    result[0] = -elt[0];
    result[1] = -elt[1];
    result[2] = -elt[2];
    result[3] = -elt[3];

    return(result);
}

inline Quat Quat::operator * (const Quat &a) const
{
    Quat result;

	Vec3 qv(elt[1],elt[2],elt[3]); Real qs = elt[0];
	Vec3 av(a[1],a[2],a[3]); Real as = a[0];

	Vec3 rv = qs*av + as*qv + cross(qv,av);
	Real rs = qs*as - dot(qv,av);

	result[1] = rv[0];
	result[2] = rv[1];
	result[3] = rv[2];
	result[0] = rs;

    return(result);
}

inline Quat Quat::operator * (Real s) const
{
    Quat result;

    result[0] = elt[0] * s;
    result[1] = elt[1] * s;
    result[2] = elt[2] * s;
    result[3] = elt[3] * s;

    return(result);
}

inline Quat Quat::operator / (Real s) const
{
    Quat result;

    result[0] = elt[0] / s;
    result[1] = elt[1] / s;
    result[2] = elt[2] / s;
    result[3] = elt[3] / s;

    return(result);
}

inline Quat operator * (Real s, const Quat &v)
{
    return(v * s);
}

// for convenience. Quat has no dot operation.
inline Real dot(const Quat &a, const Quat &b)
{
    return(a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]);
}

inline Real len(const Quat &v)
{
    return(sqrt(dot(v, v)));
}

inline Real sqrlen(const Quat &v)
{
    return(dot(v, v));
}

inline Quat norm(const Quat &v)
{
    Assert(sqrlen(v) > 0.0, "normalising length-zero vector");
    return(v / len(v));
}

inline Void normalise(Quat &v)
{
    v /= len(v);
}

inline Quat slerp (const Quat& q1, const Quat& q2, Real t)
{
	Quat result;
	Quat qq = q1;

	if (dot(qq,q2) < 0)
		qq = -q1;

	Real phi = acos(dot (qq, q2));
	Real denom = sin(phi);
	
	result = sin(phi*(1-t))/denom * qq + sin(phi*t)/denom * q2;

	return result;
}

inline Quat conjugate(const Quat &q)
{
	return Quat (q[0], -q[1], -q[2], -q[3]);
}

#endif
