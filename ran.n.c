/* Number Generator as per Numerical Recipes */

/*
This is now modified, because we want to be able to see these variables
in the main program.
*/

int iff, ma[56], inext, inextp, mj, mz, mbig;

double ran3(void) {
    static int mk, i, ii, k, mseed;


    if (iff == 0) {
        iff = 1;
        mbig = 1000000000;
        mseed = 161803398;
        mz = 0;
        mj = mseed;
        mj = mj % mbig;
        ma[55] = mj;
        mk = 1;
        for (i = 1; i < 55; i++) {
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < mz) mk = mk + mbig;
            mj = ma[ii];
        }
        for (k = 1; k <= 4; k++) {
            for (i = 1; i <= 55; i++) {
                ma[i] = ma[i] - ma[1 + ((i + 30) % 55)];
                if (ma[i] < mz) ma[i] = ma[i] + mbig;
            }
        }
        inext = 0;
        inextp = 31;
    }
    if (++inext == 56) inext = 1;
    if (++inextp == 56) inextp = 1;
    mj = ma[inext] - ma[inextp];
    if (mj < mz) mj = mj + mbig;
    ma[inext] = mj;
    return (double) mj / mbig;
}
    




