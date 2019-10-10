#include "Gillespie.hh"
#include "Simulation.hh"
#include "numtools/numtools.h"
#include <fstream>
#include <iostream>

Simulation start(int *, int *, int *, int, char *[]);

void finish();

void check(int x, int y) {
    if (x != y) {
        throw std::runtime_error("Could not read input.");
    }
}

int main(int argc, char *argv[]) {
    int nBlkEq, nBlkRun, nSteps;

    Simulation simulation = start(&nBlkEq, &nBlkRun, &nSteps, argc, argv);

    Trajectory trajectory(simulation.getNumComponents());
    simulation.run(nBlkEq, nSteps, trajectory);

    bool providedOutputFilename = false;
    for (int arg = 0; arg < argc; arg++) {
        std::string argument(argv[arg]);
        if (argument == "-o") {
            providedOutputFilename = true;
            if (arg + 1 >= argc) {
                std::cerr << "Usage: Gillespie INPUT -o OUTPUT\n\n"
                             "No filename specified after `-o` flag.\n";
                throw std::runtime_error("Invalid argument");
            }
            std::string filename(argv[arg + 1]);
            std::ofstream stream(filename);
            simulation.printTrajectory(stream, trajectory);
            arg++;
        }
    }

    if (!providedOutputFilename) {
        // we pick a default filename if none is provided on the command line
        std::ofstream stream(std::string(simulation.sys.name) + ".traj");
        simulation.printTrajectory(stream, trajectory);
    }

    finish();

    return 0;
}

Simulation start(int *nBlkEq, int *nBlkRun, int *nSteps, int argc, char *argv[]) {
    FILE *fp;
    System sys{};

    if (argc < 2) {
        std::cerr << "Usage: Gillespie INPUT -o OUTPUT\n\n"
                     "Needs at least one argument to provide input.\n";
        abort();
    }

    unsigned seed = std::mt19937::default_seed;
    auto filename = std::string();

    auto overwritten = std::vector<std::pair<std::string, std::string>>();
    auto trajectories = std::vector<std::string>();

    for (int arg = 1; arg < argc; arg++) {
        std::string argument(argv[arg]);
        if (argument == "-s") {
            if (arg + 1 >= argc) {
                std::cerr << "Usage: Gillespie INPUT -o OUTPUT -s SEED\n\n"
                             "No seed specified after `-s` flag.\n";
                abort();
            }
            std::string seedStr(argv[arg + 1]);
            seed = std::stoi(seedStr);
            arg++;
        } else if (argument == "-o") {
            // handle that later
            arg++;
        } else if (argument == "--overwrite") {
            arg++;
            if (arg >= argc) {
                std::cerr << "Usage: Gillespie INPUT --overwrite X=Y\n\n"
                             "No component specified after `--overwrite` flag.\n";
                abort();
            }
            argument = std::string(argv[arg]);
            auto index = argument.find('=');
            auto componentName = argument.substr(0, index);
            auto componentValue = argument.substr(index + 1);

            overwritten.emplace_back(componentName, componentValue);
        } else if (argument == "-t") {
            arg++;
            if (arg >= argc) {
                std::cerr << "No trajectory specified after `-t`\n";
                abort();
            }
            argument = std::string(argv[arg]);
            trajectories.push_back(argument);
        } else if (argument == "-h") {
            std::cerr << "Usage: Gillespie [options] INPUT\n\n";
            std::cerr << "Options:\n"
                         "\t-o FILE\t\tSet the output trajectory file\n"
                         "\t--overwrite X=Y\tOverwrite the initial value of component X with Y\n"
                         "\t-s SEED\t\tSpecify which seed to use\n"
                         "\t-t TRAJECTORY\tLoad a trajectory to use for the simulation\n";
            std::exit(0);
        } else {
            filename = argv[arg];
        }
    }

    if (filename.empty()) {
        fp = stdin;
        cerr << "Reading from STDIN" << std::endl;
    } else if ((fp = fopen(filename.c_str(), "r")) == NULL) {
        printf("Cannot open %s.\n", filename.c_str());
        abort();
    } else {
        cerr << "Reading input " << filename << std::endl;
    }

    char name[4096] = {0};

    check(fscanf(fp, "%4095s%*[^\n]", name), 1);
    check(fscanf(fp, "%d%*[^\n]", &sys.numComponents), 1);
    check(fscanf(fp, "%d%*[^\n]", &sys.numReactions), 1);
    check(fscanf(fp, "%d\t%d\t%d%*[^\n]", nBlkEq, nBlkRun, nSteps), 3);
    check(fscanf(fp, "%d%*[^\n]", &sys.ana), 1);

    sys.name = name;
    fclose(fp);

    printf("===============================================================================\n");
    printf("This program propagates the chemical master equation according\n");
    printf("to the Gillespie-algorithm.\n");
    printf("Log book information.\n\n");
    if (log_init(sys.name.c_str(), TRUE) != 0)
        printf("No log book information could be obtained.\n");
    printf("-------------------------------------------------------------------------------\n");
    printf("System parameters.\n\n");
    printf("Name of the run                   %8s\n", sys.name.c_str());
    printf("Number of components              %8d\n", sys.numComponents);
    printf("Number of reaction channels       %8d\n", sys.numReactions);
    printf("Number of equilibrium blocks      %8d\n", *nBlkEq);
    printf("Number of production  blocks      %8d\n", *nBlkRun);
    printf("Number of steps per block         %8d\n", *nSteps);
    printf("Frequency of analysis             %8d\n", sys.ana);

    Simulation simulation(sys, seed);
    simulation.readComponents(overwritten);
    simulation.readReactions();
    simulation.printReactions();

    for (auto &filename : trajectories) {
        ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Could not open file '" << filename << "'." << std::endl;
            std::terminate();
        }
        Trajectory trajectory;
        if ((file >> trajectory).fail()) {
            std::cerr << "Could not read in trajectory with name '" << filename << "'." << std::endl;
            std::terminate();
        }
        simulation.associateTrajectory(trajectory);
    }

    printf("===============================================================================\n");

    return simulation;
}

void finish() {
    printf("\n\n===============================================================================\n");
    printf("Run completed.\nLog book information.\n");
    log_exit();
}
