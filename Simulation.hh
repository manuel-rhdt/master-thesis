//
// Created by reinhardt on 12-09-19.
//

#ifndef GILLESPIE_SIMULATION_HH
#define GILLESPIE_SIMULATION_HH

#include "Gillespie.hh"
#include <random>

class Simulation {
private:
    std::mt19937 rnGenerator;
    std::uniform_real_distribution<double> distribution;

    System sys;

    std::vector<Reaction> reactions;
    std::vector<double> reactionRates;

    std::vector<double> componentCounts;
    std::vector<std::string> componentNames;

    void preAllocate();

    void determinePropensityFunctions(double *sumA);

    int selectReaction(double sumA);

    void updateConcentrations(int j);

public:
    explicit Simulation(System sys);

    void run(int numBlocks, int numSteps);

    void readComponents();

    void readReactions();

    void printReactions();
};

#endif //GILLESPIE_SIMULATION_HH
