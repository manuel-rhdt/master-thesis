//
// Created by reinhardt on 12-09-19.
//

#ifndef GILLESPIE_SIMULATION_HH
#define GILLESPIE_SIMULATION_HH


#include <random>
#include "Gillespie.hh"

class Simulation {
private:
    std::mt19937 rn_generator;
    std::uniform_real_distribution<double> distribution;

    System sys;
    std::vector<Reaction> reactions;
    std::vector<double> reaction_rates;
    std::vector<double> component_counts;
    std::vector<std::string> component_names;

    void pre_allocate();

public:
    Simulation(System sys);

    void run(int num_blocks, int num_steps);

    void determine_propensity_functions(double *sum_a);

    void read_components();

    void read_reactions();

    int select_reaction(double sum_a);

    void update_concentrations(int j);

    void print_reactions();
};


#endif //GILLESPIE_SIMULATION_HH
