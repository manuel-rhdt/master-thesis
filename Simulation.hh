//
// Created by reinhardt on 12-09-19.
//

#ifndef GILLESPIE_SIMULATION_HH
#define GILLESPIE_SIMULATION_HH

#include "Gillespie.hh"
#include "Trajectory.hh"
#include <random>
#include <memory>

class Simulation {
private:
    std::mt19937 rnGenerator;
    std::uniform_real_distribution<double> distribution;

    System sys;

    std::vector<Reaction> reactions;
    std::vector<double> reactionRates;
    std::vector<double> componentCounts;
    std::shared_ptr<std::vector<std::string>> componentNames;

    void preAllocate();

    void determinePropensityFunctions(double *sumA);

    int selectReaction(double sumA);

    void updateConcentrations(int j);

public:
    explicit Simulation(System sys);

    void run(int numBlocks, int numSteps, Trajectory &trajectory);

    void readComponents();

    void readReactions();

    void printReactions();

    unsigned long getNumComponents();
};

#endif //GILLESPIE_SIMULATION_HH
