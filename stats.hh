/*-----------------Header file for Gillespie simulation-------------------*/

#ifndef _STATS_H_
#define _STATS_H_

/*-----------------------Externally declared variables----------------------*/

extern statistics *Xblk, *Xrun;

/*-----------------------Externally declared functions----------------------*/

extern void blkzero(void);

extern void blkacc(double);

extern void blkout(int);

extern void statsout(void);

extern void runzero(void);

extern void runout(int);

#endif
