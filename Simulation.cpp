//
// Created by reinhardt on 12-09-19.
//

#include "Simulation.hh"
#include "Run.hh"

#include <cassert>
#include <fstream>
#include <sstream>

/**
 * Determines the stochastic time step by which to advance the simulation according to the total event rate.
 *
 * @param[in] sumA Total reaction rate for all possible reactions. Must be greater than zero or else the function will abort.
 * @param[in] randomUniformValue A random value between 0.0 and 1.0 used to determine the stochastic time step
 * @returns the determined time step
 */
double propagateTime(double sumA, double randomUniformValue)
{
    assert(sumA >= 0.0);
    if (sumA == 0.0) {
        throw std::string("Reached steady state");
    }
    return log(1.0 / randomUniformValue) / sumA;
}

/**
 * Creates a `Simulation` from a system.
 *
 * @param sys The system to initialize the simulation with.
 */
Simulation::Simulation(System sys)
    : sys(sys)
{
    rnGenerator = std::mt19937(std::mt19937::default_seed);
    distribution = std::uniform_real_distribution<double>(0.0, 1.0);
    preAllocate();
}

/**
 * Preallocate reaction and component storages.
 *
 * The reaction_rates will get initialized with a value of zero.
 */
void Simulation::preAllocate()
{
    reactions.reserve(sys.num_reactions);
    reactionRates = std::vector<double>(sys.num_reactions, 0.0);
    componentCounts.reserve(sys.num_components);
    componentNames.reserve(sys.num_components);

    //Xblk = (statistics *) calloc(sys.num_components, sizeof(statistics));
    //Xrun = (statistics *) calloc(sys.num_components, sizeof(statistics));
}

/**
 * Read components from a file called `sys.name`.components
 */
void Simulation::readComponents()
{
    std::ifstream file(std::string(sys.name) + ".components");

    std::string line;
    int lineNum = 0;
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        if (lineNum == 0) {
            int numComponents;
            lineStream >> numComponents;
            assert(numComponents == sys.num_components);
        } else {
            int count;
            std::string name;
            lineStream >> count >> name;
            componentCounts.push_back(count);
            componentNames.push_back(name);
        }
        lineNum++;
    }
    file.close();
}

/**
 * Read reactions from a file called `sys.name`.reactions
 */
void Simulation::readReactions()
{
    std::ifstream file(std::string(sys.name) + ".reactions");
    if (!file.is_open()) {
        abort();
    }

    std::string line;
    std::getline(file, line);
    std::istringstream lineStream(line);

    int numReactions;
    lineStream >> numReactions;
    assert(numReactions == sys.num_reactions);

    for (int i = 0; i < sys.num_reactions;) {
        assert(std::getline(file, line));
        std::istringstream stream1(line);
        double k; //< reaction constant
        int numReactants, numProducts;
        if (!(stream1 >> k >> numReactants >> numProducts)) {
            continue;
        }

        Reaction reaction{
            numReactants, numProducts, k
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
        i++;
    }

    file.close();
}

/**
 * Run the simulation with a given number of blocks and steps
 *
 * @param numBlocks
 * @param numSteps
 */
void Simulation::run(int runType, int numBlocks, int numSteps)
{
    Run run1(runType, &componentNames);
    for (int b = 0; b < numBlocks; b++) {
        for (int s = 0; s < numSteps; s++) {
            double sumA, rv = distribution(rnGenerator);

            determinePropensityFunctions(&sumA);
            double dt = propagateTime(sumA, rv);
            int j = selectReaction(sumA);
            updateConcentrations(j);
        }
    }
}

/**
 * Calculates the total rate.
 *
 * @param [out] sumA The total rate.
 */
void Simulation::determinePropensityFunctions(double *sumA)
{
    *sumA = 0.;
    for (int i = 0; i < sys.num_reactions; i++) {
        if (reactions[i].num_reactants == 0)
            reactionRates[i] = reactions[i].k;
        else if (reactions[i].num_reactants == 1)
            reactionRates[i] = reactions[i].k * componentCounts[reactions[i].reactant[0].index];
        else if (reactions[i].reactant[0].index == reactions[i].reactant[1].index)
            reactionRates[i] = reactions[i].k * componentCounts[reactions[i].reactant[0].index] * (componentCounts[reactions[i].reactant[1].index] - 1);
        else
            reactionRates[i] = reactions[i].k * componentCounts[reactions[i].reactant[0].index] * componentCounts[reactions[i].reactant[1].index];

        *sumA += reactionRates[i];
    }
}

/**
 * Stochastically select a reaction to perform.
 *
 * @param sumA the sum of reaction rates for all possible reactions
 * @return The index of the selected reaction.
 */
int Simulation::selectReaction(double sumA)
{
    double rs, cumulativeRate;

    rs = distribution(rnGenerator) * sumA;
    int j = 0;
    cumulativeRate = reactionRates[j];
    while (cumulativeRate < rs) {
        j++;
        cumulativeRate += reactionRates[j];
    }
    return j;
}

/**
 * Performs updates to the concentration according to the selected reaction.
 *
 * @param j The index of the selected reaction.
 */
void Simulation::updateConcentrations(int j)
{
    int i;

    for (i = 0; i < reactions[j].num_reactants; i++)
        componentCounts[reactions[j].reactant[i].index]--;

    for (i = 0; i < reactions[j].num_products; i++)
        componentCounts[reactions[j].product[i].index] += reactions[j].product[i].change;
}

void Simulation::printReactions()
{
    int i, j;

    printf("\nThe following reactions are simulated:\n\n");
    for (i = 0; i < sys.num_reactions; i++) {
        if (reactions[i].num_reactants == 0)
            printf("0");
        else
            printf("%s ", componentNames[reactions[i].reactant[0].index].c_str());
        for (j = 1; j < reactions[i].num_reactants; j++)
            printf("+ %s ", componentNames[reactions[i].reactant[j].index].c_str());
        printf(" ->  ");
        if (reactions[i].num_products == 0)
            printf("0 ");
        else
            printf("%2d %s ", reactions[i].product[0].change, componentNames[reactions[i].product[0].index].c_str());
        for (j = 1; j < reactions[i].num_products; j++)
            printf("+ %2d %s ", reactions[i].product[j].change, componentNames[reactions[i].product[j].index].c_str());
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
//            determinePropensityFunctions(&sum_a);
//            propagateTime(sum_a, &dt);
//            blkacc(dt);
//            selectReaction(&j, sum_a);
//            updateConcentrations(j);
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