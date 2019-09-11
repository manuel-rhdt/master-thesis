/*-----------------Header file for Gillespie simulation-------------------*/

#ifndef _GIL_H_
#define _GIL_H_

/*---------------------------INCLUDE'S---------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*----------------------Some uesful definitions------------------------------*/

#define EQUIL       0
#define RUN         1

#define MAXPROD     10

#define TRUE        1
#define FALSE       0
#define CANCELLED  -1
#define F_COMP(a, b) (fabs ( a-b ) < 1e-10)
#define MIN(a, b)    ( a < b ? a : b )
#define MAX(a, b)    ( a > b ? a : b )
#define EPSILON     1e-6
#define PI          3.141592653589793

/*--------------------Algebraic definitions----------------------------------*/

#define v1_v2(a, b, c) {c.x=a.x-b.x; c.y=a.y-b.y; c.z=a.z-b.z;}
#define v1v2(a, b) (a.x*b.x + a.y*b.y + a.z*b.z)
#define v1xv2(a, b, c) {c.x=a.y*b.z-a.z*b.y;c.y=a.z*b.x-a.x*b.z;c.z=a.x*b.y-a.y*b.x;}
#define cmpadd(a, b, c) {c.re = a.re+b.re; c.im = a.im + b.im;}
#define cmpxcmp(a, b, c) {c.re=a.re*b.re-a.im*b.im; c.im=a.re*b.im+a.im*b.re;}
#define cmpxcmpcon(a, b, c) {c.re=a.re*b.re+a.im*b.im; c.im=a.im*b.re-a.re*b.im;}
#define floatxcmp(a, b, c) {c.re = b.re * a; c.im = b.im * a;}
#define conjg(a, b) {b.re = a.re; b.im = -a.im;} /* b = a* */

/*-----------------------Structure definitions-------------------------------*/

typedef struct stats_type {

    int acc;

    double sum, sumsq, err, noise;

} Stats;

typedef struct sys_type {

    int Ncomp, Nreact, ana;

    double tau_blk, tau_run;

    char *name;

} Sys;

typedef struct stoch_type {

    int index, change;

} Stoch;

typedef struct reaction_type {

    int Nreact, Nprod;
    double k;
    Stoch react[2], prod[MAXPROD];

} React;

/*-----------------------Externally declared variables----------------------*/

extern Sys sys;
extern int *X;
extern char **Xname;
extern double *a;
extern React *R;

/*-----------------------Externally declared functions----------------------*/

extern void read_components(void);

extern void read_reactions(void);

extern void propagate(int);

extern void analyse(int, int, int, int);

#endif     
