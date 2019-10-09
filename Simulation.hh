//
// Created by reinhardt on 12-09-19.
//

#ifndef GILLESPIE_SIMULATION_HH
#define GILLESPIE_SIMULATION_HH

#include "Gillespie.hh"
#include "Trajectory.hh"
#include "Run.hh"
#include "Reaction.hh"
#include <random>
#include <memory>

class Simulation {
private:
    std::mt19937 rnGenerator;
    std::uniform_real_distribution<double> distribution;

    std::vector<Reaction> reactions;
    std::vector<double> propensities;
    std::vector<double> componentCounts;
    std::shared_ptr<std::vector<std::string>> componentNames;

    // TODO: document this
    // this is needed to calculate the likelihood p(x|s)
    std::vector<double> randomVariates;

    double timeStamp;

    // Some components may have associated trajectories. This means that such a component wont participate in the
    // simulation in the usual way but rather will always act as having the concentration given by the trajectory.
    //
    // It must be ensured that the entry in componentCounts for such a component is always zero.
    //
    // TODO: For now we only support a single component with associated trajectory

    std::vector<unsigned long> componentIndicesWithTrajectories;
    std::vector<Trajectory> associatedTrajectories;

    unsigned long trajectoryProgress;

    void preAllocate();

    void determinePropensityFunctions(double *sumA);

    int selectReaction();

    void updateConcentrations(int j);

public:
    explicit Simulation(System sys, unsigned seed);

    void run(int numBlocks, int numSteps, Trajectory &trajectory);


    void readReactions();

    void printReactions();

    unsigned long getNumComponents();


    bool componentHasTrajectory(unsigned long component);

    Trajectory &componentGetTrajectory(unsigned long component);

    double propagateTime();

    System sys;

    void printTrajectory(ostream &os, Trajectory &trajectory);

    void readComponents(const vector<std::pair<std::string, std::string>> &overwrite);
};

#endif //GILLESPIE_SIMULATION_HH
