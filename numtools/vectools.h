#ifndef _VECTOOLS_H_
#define _VECTOOLS_H_

/*** TYPEDEF'S ***/

typedef struct _intVector intVec;
typedef struct _dblVector dblVec;

/*** STRUCTURE'S ***/

struct _intVector {
    int x;
    int y;
    int z;
};

struct _dblVector {
    double x;
    double y;
    double z;
};

/*** PROTOTYPE'S ***/

intVec set_int_vector(int x, int y, int z);

intVec sub_int_vector(intVec a, intVec b);

intVec add_int_vector(intVec a, intVec b);

intVec rot_int_vector(intVec r, int R[3][3]);

dblVec set_dbl_vector(double x, double y, double z);

dblVec sub_dbl_vector(dblVec a, dblVec b);

dblVec add_dbl_vector(dblVec a, dblVec b);

dblVec rot_dbl_vector(dblVec r, double R[3][3]);

#endif
