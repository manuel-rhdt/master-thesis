//
// Created by reinhardt on 12-09-19.
//

#include "Simulation.hh"

#include <cassert>
#include <fstream>
#include <sstream>

/**
 * Determines the stochastic time step by which to advance the simulation according to the total event rate.
 *
 * @param[in] sum_a Total reaction rate for all possible reactions. Must be greater than zero or else the function will abort.
 * @param[in] random_uniform_value A random value between 0.0 and 1.0 used to determine the stochastic time step
 * @returns the determined time step
 */
double propagate_time(double sum_a, double random_uniform_value) {
    assert(sum_a > 0.0);
    return log(1.0 / random_uniform_value) / sum_a;
}

/**
 * Creates a `Simulation` from a system.
 *
 * @param sys The system to initialize the simulation with.
 */
Simulation::Simulation(System sys) : sys(sys) {
    rn_generator = std::mt19937(std::mt19937::default_seed);
    distribution = std::uniform_real_distribution<double>(0.0, 1.0);
    pre_allocate();
}

/**
 * Preallocate reaction and component storages.
 */
void Simulation::pre_allocate() {
    reactions.reserve(sys.num_reactions);
    reaction_rates.reserve(sys.num_reactions);
    component_counts.reserve(sys.num_components);
    component_names.reserve(sys.num_components);

    //Xblk = (Stats *) calloc(sys.num_components, sizeof(Stats));
    //Xrun = (Stats *) calloc(sys.num_components, sizeof(Stats));
}

/**
 * Read components from a file called `sys.name`.components
 */
void Simulation::read_components() {
    std::ifstream file(std::string(sys.name) + ".components");

    std::string line;
    int line_num = 0;
    while (std::getline(file, line)) {
        std::istringstream line_stream(line);
        if (line_num == 0) {
            int num_components;
            line_stream >> num_components;
            assert(num_components == sys.num_components);
        } else {
            int count;
            std::string name;
            line_stream >> count >> name;
            component_counts.push_back(count);
            component_names.push_back(name);
        }
        line_num++;
    }
    file.close();
}

/**
 * Read reactions from a file called `sys.name`.reactions
 */
void Simulation::read_reactions() {
    std::ifstream file(std::string(sys.name) + ".reactions");

    std::string line;
    std::getline(file, line);
    std::istringstream line_stream(line);

    int num_reactions;
    line_stream >> num_reactions;
    assert(num_reactions == sys.num_reactions);

    for (int i = 0; i < sys.num_reactions; i++) {
        assert(std::getline(file, line));
        std::istringstream stream1(line);
        double k; //< reaction constant
        int num_reactants, num_products;
        if (!(stream1 >> k >> num_reactants >> num_products)) {
            continue;
        };

        Reaction reaction{
                num_reactants, num_products, k
        };

        assert(std::getline(file, line));
        std::istringstream stream2(line);

        std::string dummy;
        if (reaction.num_reactants == 0) {
            stream2 >> dummy;
        } else {
            stream2 >> dummy >> reaction.reactant[0].index;
        }
        for (int j = 1; j < reaction.num_reactants; j++) {
            stream2 >> dummy >> dummy >> reaction.reactant[j].index;
        }

        stream2 >> dummy;

        if (reaction.num_products == 0) {
            stream2 >> dummy;
        } else {
            stream2 >> reaction.product[0].change >> dummy >> reaction.product[0].index;
        }
        for (int j = 1; j < reaction.num_reactants; j++) {
            stream2 >> dummy >> reaction.product[j].change >> dummy >> reaction.product[j].index;
        }

        reactions.push_back(reaction);
    }

    file.close();
}

/**
 * Run the simulation with a given number of blocks and steps
 *
 * @param num_blocks
 * @param num_steps
 */
void Simulation::run(int num_blocks, int num_steps) {
    for (int b = 0; b < num_blocks; b++) {
        for (int s = 0; s < num_steps; s++) {
            double sum_a, rv = distribution(rn_generator);

            determine_propensity_functions(&sum_a);
            double dt = propagate_time(sum_a, rv);
            int j = select_reaction(sum_a);
            update_concentrations(j);
        }
    }
}

/**
 * Calculates the total rate.
 *
 * @param [out] sum_a The total rate.
 */
void Simulation::determine_propensity_functions(double *sum_a) {
    *sum_a = 0.;
    for (int i = 0; i < sys.num_reactions; i++) {
        if (reactions[i].num_reactants == 0)
            reaction_rates[i] = reactions[i].k;
        else if (reactions[i].num_reactants == 1)
            reaction_rates[i] = reactions[i].k * component_counts[reactions[i].reactant[0].index];
        else if (reactions[i].reactant[0].index == reactions[i].reactant[1].index)
            reaction_rates[i] = reactions[i].k * component_counts[reactions[i].reactant[0].index] *
                                (component_counts[reactions[i].reactant[1].index] - 1);
        else
            reaction_rates[i] = reactions[i].k * component_counts[reactions[i].reactant[0].index] *
                                component_counts[reactions[i].reactant[1].index];

        *sum_a += reaction_rates[i];
    }
}

/**
 *
 *
 * @param sum_a the sum of reaction rates for all possible reactions
 * @return The index of the selected reaction.
 */
int Simulation::select_reaction(double sum_a) {
    double rs, cumulative_rate;

    rs = distribution(rn_generator) * sum_a;
    int j = 0;
    cumulative_rate = reaction_rates[j];
    while (cumulative_rate < rs) {
        j++;
        cumulative_rate += reaction_rates[j];
    }
    return j;
}

/**
 * Performs updates to the concentration according to the selected reaction.
 *
 * @param j The index of the selected reaction.
 */
void Simulation::update_concentrations(int j) {
    int i;

    for (i = 0; i < reactions[j].num_reactants; i++)
        component_counts[reactions[j].reactant[i].index]--;

    for (i = 0; i < reactions[j].num_products; i++)
        component_counts[reactions[j].product[i].index] += reactions[j].product[i].change;

}

void Simulation::print_reactions() {
    int i, j;

    printf("\nThe following reactions are simulated:\n\n");
    for (i = 0; i < sys.num_reactions; i++) {
        if (reactions[i].num_reactants == 0)
            printf("0");
        else
            printf("%s ", component_names[reactions[i].reactant[0].index].c_str());
        for (j = 1; j < reactions[i].num_reactants; j++)
            printf("+ %s ", component_names[reactions[i].reactant[j].index].c_str());
        printf(" ->  ");
        if (reactions[i].num_products == 0)
            printf("0 ");
        else
            printf("%2d %s ", reactions[i].product[0].change, component_names[reactions[i].product[0].index].c_str());
        for (j = 1; j < reactions[i].num_products; j++)
            printf("+ %2d %s ", reactions[i].product[j].change, component_names[reactions[i].product[j].index].c_str());
        printf("k = %4.3f\n", reactions[i].k);
    }
}

//void run(int run, int n_blk, int n_steps) {
//    int b, s, j;
//    double sum_a, dt;
//
//    runzero();
//
//    analyse(run, 0, n_steps, 0);
//
//    for (b = 0; b < n_blk; b++) {
//
//        blkzero();
//
//        for (s = 0; s < n_steps; s++) {
//            determine_propensity_functions(&sum_a);
//            propagate_time(sum_a, &dt);
//            blkacc(dt);
//            select_reaction(&j, sum_a);
//            update_concentrations(j);
//            if ((s + 1) % sys.ana == 0)
//                analyse(run, b, n_steps, s + 1);
//        }
//
//        blkout(b);
//
//        statsout();
//    }
//
//    runout(run);
//
//    return;
//}