#include "Gillespie.hh"
#include "Simulation.hh"
#include "numtools/numtools.h"
#include <fstream>

Simulation start(int *, int *, int *, int, char *[]);

void finish();

int main(int argc, char *argv[]) {
    int nBlkEq, nBlkRun, nSteps;

    Simulation simulation = start(&nBlkEq, &nBlkRun, &nSteps, argc, argv);

    Trajectory trajectory(simulation.getNumComponents());
    simulation.run(nBlkEq, nSteps, trajectory);

    std::ofstream stream("trajectory.txt");
    stream << trajectory;

    finish();

    return 0;
}

Simulation start(int *nBlkEq, int *nBlkRun, int *nSteps, int argc, char *argv[]) {
    FILE *fp;
    System sys;

    if (argc != 2) {
        printf("Needs one argument to provide input.\n");
        abort();
    }

    if ((fp = fopen(argv[1], "r")) == NULL) {
        printf("Cannot open %s.inp.\n", argv[1]);
        abort();
    }
    sys.name = (char *) calloc(30, sizeof(char));
    fscanf(fp, "%s%*s", sys.name);
    fscanf(fp, "%d%*s", &sys.numComponents);
    fscanf(fp, "%d%*s", &sys.numReactions);
    fscanf(fp, "%d\t%d\t%d%*s", nBlkEq, nBlkRun, nSteps);
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
    printf("Number of components              %8d\n", sys.numComponents);
    printf("Number of reaction channels       %8d\n", sys.numReactions);
    printf("Number of equilibrium blocks      %8d\n", *nBlkEq);
    printf("Number of production  blocks      %8d\n", *nBlkRun);
    printf("Number of steps per block         %8d\n", *nSteps);
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
