#include "Gillespie.h"
#include "stats.h"
#include "numtools/numtools.h"


/*------------------------Globally defined functions------------------------*/

void blkzero(void);

void blkacc(double);

void blkout(int);

void statsout(void);

void runzero(void);

void runout(int);

/*------------------------Globally defined variables------------------------*/

Stats *Xblk, *Xrun;


void blkzero() {
    int i;

    sys.tau_blk = 0.;
    for (i = 0; i < sys.Ncomp; i++) {
        Xblk[i].sum = 0.;
        Xblk[i].sumsq = 0.;
        Xblk[i].acc = 0;
    }

    return;
}

void blkacc(double dt) {
    int i;

    sys.tau_blk += dt;
    for (i = 0; i < sys.Ncomp; i++) {
        Xblk[i].sum += X[i] * dt;
        Xblk[i].sumsq += X[i] * X[i] * dt;
        Xblk[i].acc++;
    }

    return;
}

void blkout(int block) {
    int i;

    for (i = 0; i < sys.Ncomp; i++) {
        Xblk[i].sum /= sys.tau_blk;
        Xblk[i].sumsq /= sys.tau_blk;
        Xblk[i].err = Xblk[i].sumsq - Xblk[i].sum * Xblk[i].sum; /*sigma^2*/
        if (Xblk[i].err > 0.) Xblk[i].err = sqrt(Xblk[i].err); /*sigma*/
        if (Xblk[i].sum != 0.)
            Xblk[i].noise = Xblk[i].err / Xblk[i].sum;
        else
            Xblk[i].noise = 0.;
    }

    sys.tau_run += sys.tau_blk;
    for (i = 0; i < sys.Ncomp; i++) {
        Xrun[i].sum += Xblk[i].sum;
        Xrun[i].sumsq += Xblk[i].sum * Xblk[i].sum;
        Xrun[i].noise += Xblk[i].noise;
        Xrun[i].acc++;
    }


    printf("\n\nThe results of block %d\n\n", block);
    printf("Elapsed time         %8.4f\n\n", sys.tau_blk);

    printf("Component\tblock average\tdeviation\t   noise\t #counts\n");
    for (i = 0; i < sys.Ncomp; i++) {
        printf("%s\t\t%8.4f\t%8.4f\t%8.4f\t%8d\n",
               Xname[i],
               Xblk[i].sum,
               Xblk[i].err,
               Xblk[i].noise,
               Xblk[i].acc);
    }

    return;
}

void runzero(void) {
    int i;

    sys.tau_run = 0.;
    for (i = 0; i < sys.Ncomp; i++) {
        Xrun[i].sum = 0.;
        Xrun[i].sumsq = 0.;
        Xrun[i].acc = 0;
        Xrun[i].noise = 0.;
    }

    return;
}

void statsout(void) {
    int i;
    Stats stats_tp;
    FILE *fp;
    char filename[40];

    sprintf(filename, "%s.so", sys.name);
    fp = fopen(filename, "w");
    for (i = 0; i < sys.Ncomp; i++) {
        stats_tp.sum = Xrun[i].sum / Xrun[i].acc;
        stats_tp.sumsq = Xrun[i].sumsq / Xrun[i].acc;
        stats_tp.err = (stats_tp.sumsq - stats_tp.sum * stats_tp.sum) /
                       (double) Xrun[i].acc;
        if (stats_tp.err > 0.) stats_tp.err = sqrt(stats_tp.err);
        stats_tp.noise = Xrun[i].noise / Xrun[i].acc;

        fprintf(fp, "%s\t\t%8.4f\t%8.4f\t%8.4f\t%8d\n", Xname[i], stats_tp.sum, stats_tp.err, stats_tp.noise,
                stats_tp.acc);
    }
    fclose(fp);

    return;
}

void runout(int run) {
    int i;

    printf("\n==============================================================================\n\n");
    if (run == EQUIL)
        printf("The averages of the equilibration run.\n\n");
    else
        printf("The average of the production run.\n\n");

    printf("Elapsed time         %8.4f\n\n", sys.tau_run);

    printf("Run averages.\n\n");
    printf("Component\tRun average\t   error\t   noise\t\t#counts\n");
    for (i = 0; i < sys.Ncomp; i++) {
        Xrun[i].sum /= Xrun[i].acc;
        Xrun[i].sumsq /= Xrun[i].acc;
        Xrun[i].err = (Xrun[i].sumsq - Xrun[i].sum * Xrun[i].sum) /
                      (double) Xrun[i].acc;
        if (Xrun[i].err > 0.) Xrun[i].err = sqrt(Xrun[i].err);
        Xrun[i].noise /= Xrun[i].acc;
        printf("%s\t\t%8.4f\t%8.4f\t%8.4f\t%8d\n", Xname[i], Xrun[i].sum, Xrun[i].err, Xrun[i].noise, Xrun[i].acc);
    }

    return;
}
