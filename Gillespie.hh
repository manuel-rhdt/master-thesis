/*-----------------Header file for Gillespie simulation-------------------*/

#ifndef _GIL_H_
#define _GIL_H_

/*---------------------------INCLUDE'S---------------------------------------*/

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>

/*----------------------Some useful definitions------------------------------*/

#define EQUIL       0
#define RUN         1

#define MAXPROD     10

#define TRUE        1
#define FALSE       0
#define PI          3.141592653589793

/*-----------------------Structure definitions-------------------------------*/

typedef struct statsType {
    int acc;
    double sum, sumsq, err, noise;
} statistics;

typedef struct sys_type {

    int num_components, num_reactions, ana;

    double tau_blk, tau_run;

    char *name;

} System;

typedef struct stoch_type {

    int index, change;

} Stoch;

/**
 * A type to describe a single reaction (?)
 */
typedef struct reaction_type {

    int num_reactants, num_products;
    double k;
    Stoch reactant[2], product[MAXPROD];

} Reaction;

/*-----------------------Externally declared functions----------------------*/

extern void read_components(void);

extern void read_reactions(void);

extern void analyse(int, int, int, int);

#endif     
