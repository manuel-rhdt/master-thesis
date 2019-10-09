/*-----------------Header file for Gillespie simulation-------------------*/

#ifndef _GIL_H_
#define _GIL_H_

/*---------------------------INCLUDE'S---------------------------------------*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

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

class System {
public:
    int numComponents{0}, numReactions{0}, ana{0};
    std::string name;
};

typedef struct stoch_type {

    int index, change;

} Stoch;


extern void analyse(int, int, int, int);

#endif
