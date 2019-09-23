/*-----------------Header file for Gillespie simulation-------------------*/

#ifndef _GIL_H_
#define _GIL_H_

/*---------------------------INCLUDE'S---------------------------------------*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

/*----------------------Some useful definitions------------------------------*/

#define EQUIL 0
#define RUN 1

#define MAXPROD 10

#define TRUE 1
#define FALSE 0
#define PI 3.141592653589793

/*-----------------------Structure definitions-------------------------------*/

typedef struct statsType {
    int acc;
    double sum, sumsq, err, noise;
} statistics;

typedef struct sys_type {

    int numComponents, numReactions, ana;

    char *name;

} System;

typedef struct stoch_type {

    int index, change;

} Stoch;


extern void analyse(int, int, int, int);

#endif
