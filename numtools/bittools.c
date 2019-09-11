/* bittools.c
   Tools for efficient memory handling, using bitmaps.
   Joel Wijngaarde, 1998
   */

/*** MAIN INCLUDE ***/

#include "bittools.h"

/*---------------------------------------------------------------------------*/

bitfield *malloc_3D_bitfield(unsigned int x, unsigned int y, unsigned int z) {
    bitfield *tmp = NULL;
    int bits, nb;

    nb = sizeof(unsigned long int) * 8;
    tmp = (bitfield *) malloc(sizeof(bitfield));
    if (!tmp) {
        perror("malloc_3D_bitfield");
        exit(-1);
    }
    bits = (x * y * z);
    bits = (bits / nb) + (bits % nb ? 1 : 0);

    tmp->type = BITFIELD3D;
    tmp->sizex = x;
    tmp->sizey = y;
    tmp->sizez = z;
    tmp->nbits = nb;
    if (nb == 16) /* not machine independend */
        tmp->shft = 3;
    else if (nb == 32)
        tmp->shft = 4;
    else if (nb == 64)
        tmp->shft = 5;
    tmp->mask = ~(0 << tmp->shft);
    tmp->data = (unsigned long int *) malloc(sizeof(unsigned long int) * bits);
    if (!tmp->data) {
        perror("malloc_3D_bitfield");
        exit(-1);
    }
    return tmp;
}

void free_3D_bitfield(bitfield *bf) {
    free(bf->data);
    free(bf);
}


/* This should be MACRO's, so use MACRO's instead! These functions are 
   !NOT! up to date, and use slower algoritmes */
char bit_query(bitfield *bf, unsigned const int x, unsigned const int y,
               unsigned const int z) {
    unsigned const int bit = z + (bf->sizez) * (y + (bf->sizey) * x);
    return (bf->data[bit / LONG_BIT] & MASK(bit % LONG_BIT));
}

char bit_set(bitfield *bf, unsigned const int x, unsigned const int y,
             unsigned const int z) {
    unsigned const int bit = z + (bf->sizez) * (y + (bf->sizey) * x);
    (bf->data[bit / LONG_BIT] |= MASK(bit % LONG_BIT));
    return 1;
}

char bit_clear(bitfield *bf, unsigned const int x, unsigned const int y,
               unsigned const int z) {
    unsigned const int bit = z + (bf->sizez) * (y + (bf->sizey) * x);
    (bf->data[bit / LONG_BIT] &= ~MASK(bit % LONG_BIT));
    return 1;
}
