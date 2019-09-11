#include "propagate.h"
#include "Gillespie.h"
#include "stats.h"

/*------------------------Globally defined functions-------------------------*/

void propagate(int);

/*---------------------- Externally defined functions-----------------------*/

extern double ran3(void);

extern int iff, ma[56], inext, inextp, mj, mz, mbig;

/*------------------------Locally defined functions-------------------------*/

void determine_propensity_functions(double *);

void propagate_time(double, double *);

void select_reaction(int *, double);

void update_concentrations(int);


void run(int run, int n_blk, int n_steps) {
    int b, s, j;
    double sum_a, dt;

    runzero();

    analyse(run, 0, n_steps, 0);

    for (b = 0; b < n_blk; b++) {

        blkzero();

        for (s = 0; s < n_steps; s++) {

            determine_propensity_functions(&sum_a);

            propagate_time(sum_a, &dt);

            blkacc(dt);

            select_reaction(&j, sum_a);

            update_concentrations(j);

            if ((s + 1) % sys.ana == 0) analyse(run, b, n_steps, s + 1);
        }

        blkout(b);

        statsout();

    }

    runout(run);

    return;
}

void determine_propensity_functions(double *sum_a) {
    int i;

    *sum_a = 0.;
    for (i = 0; i < sys.Nreact; i++) {

        if (R[i].Nreact == 0)
            a[i] = R[i].k;
        else if (R[i].Nreact == 1)
            a[i] = R[i].k * X[R[i].react[0].index];
        else if (R[i].react[0].index == R[i].react[1].index)
            a[i] = R[i].k * X[R[i].react[0].index] * (X[R[i].react[1].index] - 1);
        else
            a[i] = R[i].k * X[R[i].react[0].index] * X[R[i].react[1].index];

        *sum_a += a[i];
    }

    return;
}

void propagate_time(double sum_a, double *dt) {
    int i;
    if (sum_a > 0.) {
        *dt = log(1. / ran3()) / sum_a;
    } else {
        printf("Not a single reaction can occur.\n");
        printf("The run will be terminated.\n");
        printf("sum_a is %f\n", sum_a);
        for (i = 0; i < sys.Ncomp; i++) printf("%d: N is %d\n", i, X[i]);
        abort();
    }

}

void select_reaction(int *j, double sum_a) {
    double rs, cumu_a;

    rs = ran3() * sum_a;
    *j = 0;
    cumu_a = a[*j];
    while (cumu_a < rs) {
        (*j)++;
        cumu_a += a[*j];
    }
    return;
}

void update_concentrations(int j) {
    int i, k;

    for (i = 0; i < R[j].Nreact; i++)
        X[R[j].react[i].index]--;

    for (i = 0; i < R[j].Nprod; i++)
        X[R[j].prod[i].index] += R[j].prod[i].change;

    return;
}
