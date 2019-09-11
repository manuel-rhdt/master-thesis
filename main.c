#include "Gillespie.h"
#include "stats.h"
#include "numtools/numtools.h"
#include "propagate.h"

/*------------------------Globally defined variables------------------------*/

Sys sys;
int *X;
char **Xname;
double *a, sum_a;
React *R;

/*------------------------Locally defined functions--------------------------*/

void start(int *, int *, int *);

void allocate_memory();

void print_reactions();

void finish();

int main(void) {
    int n_blk_eq, n_blk_run, n_steps;

    start(&n_blk_eq, &n_blk_run, &n_steps);

    run(EQUIL, n_blk_eq, n_steps);

    run(RUN, n_blk_run, n_steps);

    finish();

    return 0;
}

void start(int *n_blk_eq, int *n_blk_run, int *n_steps) {
    int i, j;
    FILE *fp;

    if ((fp = fopen("Gillespie.inp", "r")) == NULL) {
        printf("Cannot open Gillespie.inp.\n");
        abort();
    }
    sys.name = (char *) calloc(30, sizeof(char));
    fscanf(fp, "%s%*s", sys.name);
    fscanf(fp, "%d%*s", &sys.Ncomp);
    fscanf(fp, "%d%*s", &sys.Nreact);
    fscanf(fp, "%d\t%d\t%d%*s", n_blk_eq, n_blk_run, n_steps);
    fscanf(fp, "%d%*s", &sys.ana);
    fclose(fp);

    printf("===============================================================================\n");
    printf("This program propagates the chemical master equation according\n");
    printf("to the Gillespie-algorithm.\n");
    printf("Log book information.\n\n");
    if (log_init(sys.name, TRUE) != 0)
        printf("No log book information could be obtained.\n");
    printf("-------------------------------------------------------------------------------\n");
    printf("System parameters.\n\n");
    printf("Name of the run                   %8s\n", sys.name);
    printf("Number of components              %8d\n", sys.Ncomp);
    printf("Number of reaction channels       %8d\n", sys.Nreact);
    printf("Number of equilibrium blocks      %8d\n", *n_blk_eq);
    printf("Number of production  blocks      %8d\n", *n_blk_run);
    printf("Number of steps per block         %8d\n", *n_steps);
    printf("Frequency of analysis             %8d\n", sys.ana);

    allocate_memory();

    read_components();

    read_reactions();

    print_reactions();

    printf("===============================================================================\n");
}

void allocate_memory() {
    int i;

    R = (React *) calloc(sys.Nreact, sizeof(React));
    a = (double *) calloc(sys.Nreact, sizeof(double));

    X = (int *) calloc(sys.Ncomp, sizeof(int));
    Xname = (char **) calloc(sys.Ncomp, sizeof(char *));
    for (i = 0; i < sys.Ncomp; i++) Xname[i] = (char *) calloc(30, sizeof(char));
    Xblk = (Stats *) calloc(sys.Ncomp, sizeof(Stats));
    Xrun = (Stats *) calloc(sys.Ncomp, sizeof(Stats));

    return;
}


void print_reactions() {
    int i, j;

    printf("\nThe following reactions are simulated:\n\n");
    for (i = 0; i < sys.Nreact; i++) {
        if (R[i].Nreact == 0)
            printf("0");
        else
            printf("%s ", Xname[R[i].react[0].index]);
        for (j = 1; j < R[i].Nreact; j++)
            printf("+ %s ", Xname[R[i].react[j].index]);
        printf(" ->  ");
        if (R[i].Nprod == 0)
            printf("0 ");
        else
            printf("%2d %s ", R[i].prod[0].change, Xname[R[i].prod[0].index]);
        for (j = 1; j < R[i].Nprod; j++)
            printf("+ %2d %s ", R[i].prod[j].change, Xname[R[i].prod[j].index]);
        printf("k = %4.3f\n", R[i].k);
    }
    return;
}

void finish(void) {
    printf("\n\n===============================================================================\n");
    printf("Run completed.\nLog book information.\n");
    log_exit();
}
