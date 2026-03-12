#include "Quat.h"
//#include "svl/Mat3.h"
#include <cctype>
#include <iomanip>

//
// Quat in SVL: [(q0,q1,q2),q3]
// update: 
// Modify to conform with most literature:
// Quat in SVL [q0, (q1,q2,q3)]
// (may need to modify other SVL files relating to quat!!!
//
// implement:
//
// Quaterion multiplication
// 3x3 -> mQuat; 
// mQuat->3x3 (NY) .... already available 
// slerp
//
// find application that demonstrate the needs ...
//

Mat3 Rot3 (const Quat &q)
{
	Mat3 result; 
	Vec4 v(q[0],q[1],q[2],q[3]);
 	result.MakeRot(v); 
	return(result); 
}

Mat4 HRot4 (const Quat &q)
{
	Mat4 result; 
	Vec4 v(q[0],q[1],q[2],q[3]);
 	result.MakeHRot(v);
	return(result); 
}

Quat::Quat(const Mat3 & m)
{   
	Real compare[4], max;
    compare[0] = 1 + m[0][0] + m[1][1] + m[2][2];
    compare[1] = 1 + m[0][0] - m[1][1] - m[2][2];
    compare[2] = 1 - m[0][0] + m[1][1] - m[2][2];
    compare[3] = 1 - m[0][0] - m[1][1] + m[2][2];

	int i, which;
	for (i = 0, which = -1, max = 0; i < 4; i++)
	if (fabs(compare[i]) > max)
	    which = i, max = fabs(compare[i]);

	Real q0,q1,q2,q3;
	switch (which) {
      case 0:
		q0 = 0.5*sqrt (compare[which]);
		q1 = 0.25*(m[2][1] - m[1][2])/q0;
		q2 = 0.25*(m[0][2] - m[2][0])/q0;
		q3 = 0.25*(m[1][0] - m[0][1])/q0;
		break;

      case 1:
		q1 = 0.5*sqrt (compare[which]);
		q0 = 0.25*(m[2][1] - m[1][2])/q1;
		q2 = 0.25*(m[0][1] + m[1][0])/q1;
		q3 = 0.25*(m[0][2] + m[2][0])/q1;
		break;

      case 2:
		q2 = 0.5*sqrt (compare[which]);
		q0 = 0.25*(m[0][2] - m[2][0])/q2;
		q1 = 0.25*(m[0][1] + m[1][0])/q2;
		q3 = 0.25*(m[1][2] + m[2][1])/q2;
		break;

      case 3:
		q3 = 0.5*sqrt (compare[which]);
		q0 = 0.25*(m[1][0] - m[0][1])/q3;
		q1 = 0.25*(m[0][2] + m[2][0])/q3;
		q2 = 0.25*(m[1][2] + m[2][1])/q3;
		break;
    }
	elt[1] = q1, elt[2] = q2, elt[3] = q3, elt[0] = q0;
}


ostream &operator << (ostream &s, const Quat &v)
{
    Int w = s.width();

    return(s << '[' << v[0] << ' ' << setw(w) << v[1] << ' '
        << setw(w) << v[2] << ' ' << setw(w) << v[3] << ']');
}

istream &operator >> (istream &s, Quat &v)
{
    Quat    result;
    Char    c;

    // Expected format: [1 2 3 4]

    while (s >> c && isspace(c))
        ;

    if (c == '[')
    {
        s >> result[0] >> result[1] >> result[2] >> result[3];

        if (!s)
        {
            cerr << "Error: Expected number while reading vector\n";
            return(s);
        }

        while (s >> c && isspace(c))
            ;

        if (c != ']')
        {
            s.clear(ios::failbit);
            cerr << "Error: Expected ']' while reading vector\n";
            return(s);
        }
    }
    else
    {
        s.clear(ios::failbit);
        cerr << "Error: Expected '[' while reading vector\n";
        return(s);
    }

    v = result;
    return(s);
}

