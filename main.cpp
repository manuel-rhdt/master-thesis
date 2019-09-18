#include "Gillespie.hh"
#include "Simulation.hh"
#include "numtools/numtools.h"
#include "stats.hh"
#include <random>
#include <fstream>

Simulation start(int *, int *, int *);

void finish();

int main(void) {
    int n_blk_eq, n_blk_run, n_steps;

    Simulation simulation = start(&n_blk_eq, &n_blk_run, &n_steps);

    Trajectory trajectory;
    std::ifstream inputStream("int_traj.txt");
    inputStream >> trajectory;
//    simulation.run(n_blk_eq, n_steps, trajectory);

    std::ofstream stream("int2_traj.txt");
    stream << trajectory;

    finish();

    return 0;
}

Simulation start(int *n_blk_eq, int *n_blk_run, int *n_steps) {
    FILE *fp;
    System sys;

    if ((fp = fopen("Gillespie.inp", "r")) == NULL) {
        printf("Cannot open Gillespie.inp.\n");
        abort();
    }
    sys.name = (char *) calloc(30, sizeof(char));
    fscanf(fp, "%s%*s", sys.name);
    fscanf(fp, "%d%*s", &sys.num_components);
    fscanf(fp, "%d%*s", &sys.num_reactions);
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
    printf("Number of components              %8d\n", sys.num_components);
    printf("Number of reaction channels       %8d\n", sys.num_reactions);
    printf("Number of equilibrium blocks      %8d\n", *n_blk_eq);
    printf("Number of production  blocks      %8d\n", *n_blk_run);
    printf("Number of steps per block         %8d\n", *n_steps);
    printf("Frequency of analysis             %8d\n", sys.ana);

    Simulation simulation(sys);
    simulation.readComponents();
    simulation.readReactions();
    simulation.printReactions();

    printf("===============================================================================\n");

    return simulation;
}

void finish(void) {
    printf("\n\n===============================================================================\n");
    printf("Run completed.\nLog book information.\n");
    log_exit();
}
