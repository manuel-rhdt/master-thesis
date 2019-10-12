
#include "vectools.h"

/*---------------------------------------------------------------------------*/


/* set_int_vector():
** initialize a Vector with the specified integer values for x, y and z.
*/
#ifdef __GNUC__

__inline__ intVec set_int_vector(int x, int y, int z)
#else
intVec set_int_vector( int x, int y, int z )
#endif
{
    intVec r;
    r.x = x;
    r.y = y;
    r.z = z;
    return r;
}
/* sub_int_vector():
** subtracts two integer vectors.
*/
#ifdef __GNUC__

__inline__ intVec sub_int_vector(intVec a, intVec b)
#else
intVec sub_int_vector( intVec a, intVec b )
#endif
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}
/* add_int_vector():
** adds two integer vectors.
*/
#ifdef __GNUC__

__inline__ intVec add_int_vector(intVec a, intVec b)
#else
intVec add_int_vector( intVec a, intVec b )
#endif
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}
/* rot_int_vector():
** rotate a float vector using a certain rotataion matrix 'R'
** containing ints.
*/
#ifdef __GNUC__

__inline__ intVec rot_int_vector(intVec r, int R[3][3])
#else
intVec rot_int_vector( intVec r, int R[3][3] )
#endif
{
    return set_int_vector(R[0][0] * r.x + R[0][1] * r.y + R[0][2] * r.z,
                          R[1][0] * r.x + R[1][1] * r.y + R[1][2] * r.z,
                          R[2][0] * r.x + R[2][1] * r.y + R[2][2] * r.z);
}
/* set_dbl_vector():
** initialize a Vector with the specified double values for x, y and z.
*/
#ifdef __GNUC__

__inline__ dblVec set_dbl_vector(double x, double y, double z)
#else
dblVec set_dbl_vector( double x, double y, double z )
#endif
{
    dblVec r;
    r.x = x;
    r.y = y;
    r.z = z;
    return r;
}
/* sub_dbl_vector():
** subtracts two float vectors.
*/
#ifdef __GNUC__

__inline__ dblVec sub_dbl_vector(dblVec a, dblVec b)
#else
dblVec sub_dbl_vector( dblVec a, dblVec b )
#endif
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}
/* add_dbl_Vector():
** adds two float vectors.
*/
#ifdef __GNUC__

__inline__ dblVec add_dbl_vector(dblVec a, dblVec b)
#else
dblVec add_dbl_vector( dblVec a, dblVec b )
#endif
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}
/* rot_dbl_vector():
** rotate a float vector using a certain rotataion matrix 'R'
** conatining doubles.
*/
#ifdef __GNUC__

__inline__ dblVec rot_dbl_vector(dblVec r, double R[3][3])
#else
dblVec rot_dbl_vector( dblVec r, double R[3][3] )
#endif
{
    return set_dbl_vector(R[0][0] * r.x + R[0][1] * r.y + R[0][2] * r.z,
                          R[1][0] * r.x + R[1][1] * r.y + R[1][2] * r.z,
                          R[2][0] * r.x + R[2][1] * r.y + R[2][2] * r.z);
}


