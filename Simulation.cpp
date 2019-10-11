//
// Created by reinhardt on 12-09-19.
//

#include "Simulation.hh"
#include "Run.hh"
#include "Block.hh"

#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>

/**
 * Determines the stochastic time step by which to advance the simulation according to the total event rate.
 *
 * @param[in] sumA Total reaction rate for all possible reactions. Must be greater than zero or else the function will abort.
 * @param[in] randomUniformValue A random value between 0.0 and 1.0 used to determine the stochastic time step
 * @returns the determined time step
 */
double propagateTime(double sumA, double randomUniformValue) {
    assert(sumA >= 0.0);
    if (sumA == 0.0) {
        throw std::string("Reached steady state");
    }
    return -log(randomUniformValue) / sumA;
}

/**
 * Creates a `Simulation` from a system.
 *
 * @param sys The system to initialize the simulation with.
 */
Simulation::Simulation(System sys, unsigned seed)
        : sys(sys), trajectoryProgress(0), timeStamp(0.0) {
    rnGenerator = std::mt19937(seed);
    distribution = std::uniform_real_distribution<double>(0.0, 1.0);
    reactions = std::vector<Reaction>();
    propensities = std::vector<double>();
    componentCounts = std::vector<double>();
    componentNames = std::make_shared<std::vector<std::string>>(std::vector<std::string>());
    preAllocate();
    componentIndicesWithTrajectories = std::vector<unsigned long>();
    associatedTrajectories = std::vector<Trajectory>();
    randomVariates = std::vector<double>();
}

/**
 * Preallocate reaction and component storages.
 *
 * The reaction_rates will get initialized with a value of zero.
 */
void Simulation::preAllocate() {
    reactions.reserve(sys.numReactions);
    propensities = std::vector<double>(sys.numReactions, 0.0);
    componentCounts.reserve(sys.numComponents);
    componentNames->reserve(sys.numComponents);
}

/**
 * Read components from a file called `sys.name`.components
 *
 * Allows overriding initial values for specific components from the command line through the `overwrite` parameter.
 */
void Simulation::readComponents(const std::vector<std::pair<std::string, std::string>> &overwrite) {
    std::ifstream file(std::string(sys.name) + ".components");

    std::string line;
    int lineNum = 0;
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        if (lineNum == 0) {
            int numComponents;
            lineStream >> numComponents;
            assert(numComponents == sys.numComponents);
        } else {
            // decide whether we only know the initial concentration of a component or its whole trajectory
            std::string initialValue;
            std::string name;
            if (lineStream >> initialValue >> name) {
                for (auto override : overwrite) {
                    if (override.first == name) {
                        initialValue = override.second;
                    }
                }

                if (initialValue[0] == '@') {
                    componentCounts.push_back(0.0);
                    componentNames->push_back(name);
                    std::string filename = initialValue.substr(1);

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

                    associateTrajectory(trajectory);
                } else {
                    auto count = std::stoi(initialValue);
                    componentCounts.push_back(count);
                    componentNames->push_back(name);
                }
            }
        }
        lineNum++;
    }
    file.close();
}

void Simulation::associateTrajectory(Trajectory trajectory) {
    for (long component = 0; component < componentNames->size(); component++) {
        auto &tn = *trajectory.componentNames;

        if (std::find(tn.begin(), tn.end(), (*componentNames)[component]) != tn.end()) {
            // trajectory contains a component which is also contained in the simulation
            std::cout << "\nAssociating component '" << (*componentNames)[component] << "' with a trajectory\n\n";
            componentIndicesWithTrajectories.push_back(component);
            associatedTrajectories.push_back(trajectory);
        }
    }
}

/**
 * Determines the stochastic time step by which to advance the simulation.
 *
 * @return The chosen time step.
 */
double Simulation::propagateTime() {
    double uniformVariate = distribution(rnGenerator);
    double negLogUnifVal = -log(uniformVariate);
    double accumulatedTime = 0.0;

    randomVariates.push_back(uniformVariate);

    while (true) {
        double totalPropensity = 0.0;

        // sum up the propensities for the current trajectory state
        for (int i = 0; i < sys.numReactions; i++) {
            propensities[i] = reactions[i].reactionConstant;

            // TODO: x * (x - 1) * (x - 2) * ...
            for (auto reactant : reactions[i].reactants) {
                if (componentHasTrajectory(reactant.index)) {
                    propensities[i] *= associatedTrajectories[0].componentCounts[0][trajectoryProgress];
                } else {
                    propensities[i] *= componentCounts[reactant.index];
                }
            }

            if (propensities[i] < 0.0) {
                throw std::runtime_error("Error: trajectory contains negative values.");
            }

            totalPropensity += propensities[i];
        }

        if (totalPropensity == 0.0) {
            throw std::runtime_error("System reached steady state.");
        }

        // There is a maximum time step for which `totalPropensity` is valid. After this time we need to recalculate the
        // propensities to take into account changes in the external trajectory.
        double maxTimeStep;
        if (associatedTrajectories.empty() || trajectoryProgress == associatedTrajectories[0].timeStamps.size() - 1) {
            // The propensities do not depend on any external trajectory.
            maxTimeStep = std::numeric_limits<double>::infinity();
        } else {
            if (timeStamp > associatedTrajectories[0].timeStamps[trajectoryProgress]) {
                maxTimeStep = associatedTrajectories[0].timeStamps[trajectoryProgress + 1] - timeStamp;
            } else {
                maxTimeStep = associatedTrajectories[0].timeStamps[trajectoryProgress + 1] -
                              associatedTrajectories[0].timeStamps[trajectoryProgress];
            }
        }

        if (maxTimeStep * totalPropensity < negLogUnifVal) {
            negLogUnifVal -= maxTimeStep * totalPropensity;
            // ensure that trajectoryProgress always stays in bounds
            // this should be the only time we change trajectoryProgress
            trajectoryProgress = std::min(trajectoryProgress + 1, associatedTrajectories[0].timeStamps.size() - 1);
            accumulatedTime += maxTimeStep;
            continue;
        } else {
            auto timeStep = negLogUnifVal / totalPropensity + accumulatedTime;
            return timeStep;
        }

    }


}


bool Simulation::componentHasTrajectory(unsigned long component) {
    for (auto compWithTraj : componentIndicesWithTrajectories) {
        if (component == compWithTraj) {
            return true;
        }
    }
    return false;
}

Trajectory &Simulation::componentGetTrajectory(unsigned long component) {
    for (int i = 0; i < componentIndicesWithTrajectories.size(); ++i) {
        if (component == componentIndicesWithTrajectories[i]) {
            return associatedTrajectories[i];
        }
    }
    throw std::invalid_argument("Component has no associated trajectory.");
}

/**
 * Read reactions from a file called `sys.name`.reactions
 */
void Simulation::readReactions() {
    std::ifstream file(std::string(sys.name) + ".reactions");
    if (!file.is_open()) {
        throw std::runtime_error("Could not open '" + std::string(sys.name) + ".reactions" + "'");
    }

    std::string line;
    std::getline(file, line);
    std::istringstream lineStream(line);

    int numReactions = 0;
    lineStream >> numReactions;
    assert(numReactions == sys.numReactions);

    for (int i = 0; i < sys.numReactions;) {
        if (!std::getline(file, line)) {
            std::cerr << "Couldn't read " << sys.numReactions << " reactions";
            std::terminate();
        }
        std::istringstream stream1(line);
        double k = 0.0; //< reaction constant
        int numReactants = 0, numProducts = 0;
        if (!(stream1 >> k >> numReactants >> numProducts)) {
            if (file.eof()) {
                std::terminate();
            }
            line.clear();
            continue;
        }

        Reaction reaction;
        reaction.reactionConstant = k;

        if (!std::getline(file, line)) {
            std::cerr << "Couldn't read " << sys.numReactions << " reactions";
            std::terminate();
        }
        std::istringstream stream2(line);

        std::string dummy;
        if (numReactants == 0) {
            stream2 >> dummy;
        } else {
            reaction.reactants.push_back(ReactionComponent{0, 0});
            stream2 >> dummy >> reaction.reactants[0].index;
        }
        for (int j = 1; j < numReactants; j++) {
            reaction.reactants.push_back(ReactionComponent{0, 0});
            stream2 >> dummy >> dummy >> reaction.reactants[j].index;
        }

        stream2 >> dummy;

        if (numProducts == 0) {
            stream2 >> dummy;
        } else {
            reaction.products.push_back(ReactionComponent{0, 0});
            stream2 >> reaction.products[0].change >> dummy >> reaction.products[0].index;
        }
        for (int j = 1; j < numProducts; j++) {
            reaction.products.push_back(ReactionComponent{0, 0});
            stream2 >> dummy >> reaction.products[j].change >> dummy >> reaction.products[j].index;
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
 * @param trajectory The trajectory where the simulated results will be stored. The trajectory must have been initialized
 * with the correct number of components.
 */
void Simulation::run(int numBlocks, int numSteps, Trajectory &trajectory) {
    Run run(sys.numComponents, *componentNames);
    trajectory.componentNames = componentNames;
    for (int b = 0; b < numBlocks; b++) {
        Block block(sys.numComponents, componentNames);
        for (int s = 0; s < numSteps; s++) {
            // sumA = 0.0, rv = distribution(rnGenerator);
            //determinePropensityFunctions(&sumA);
            double dt = propagateTime();
            timeStamp += dt;
            block.accumulate(dt, componentCounts);
            int j = selectReaction();
            updateConcentrations(j);

            trajectory.insertTimestamp(dt);
            trajectory.reactions.back() = j;
            for (unsigned long comp = 0; comp < sys.numComponents; comp++) {
                trajectory.componentCounts[comp].back() = componentCounts[comp];
            }
        }
        std::cout << block;
        block.updateRun(run);
    }
    run.finish();
}

/**
 * Calculates the total rate.
 *
 * @param [out] sumA The total rate.
 */
void Simulation::determinePropensityFunctions(double *sumA) {
    *sumA = 0.;
    for (int i = 0; i < sys.numReactions; i++) {
        propensities[i] = reactions[i].reactionConstant;

        // TODO: x * (x - 1) * (x - 2) * ...
        for (auto reactant : reactions[i].reactants) {
            propensities[i] *= componentCounts[reactant.index];
        }

        *sumA += propensities[i];
    }
}

/**
 * Stochastically select a reaction to perform.
 *
 * @param sumA the sum of reaction rates for all possible reactions
 * @return The index of the selected reaction.
 */
int Simulation::selectReaction() {
    double rs, cumulativeRate;

    double sumA = 0.0;
    for (auto val : propensities) {
        sumA += val;
    }

    rs = distribution(rnGenerator) * sumA;
    int j = 0;
    cumulativeRate = propensities[j];
    while (cumulativeRate < rs) {
        j++;
        cumulativeRate += propensities[j];
    }
    return j;
}

/**
 * Performs updates to the concentration according to the selected reaction.
 *
 * @param j The index of the selected reaction.
 */
void Simulation::updateConcentrations(int j) {
    int i;

    for (i = 0; i < reactions[j].reactants.size(); i++) {
        auto reactant = reactions[j].reactants[i].index;
        // If the component concentration is derived from a trajectory, we just set it to the correct value.
        if (componentHasTrajectory(reactant)) {
            auto &name = (*componentNames)[reactant];
            componentCounts[reactant] = componentGetTrajectory(reactant).getComponent(name)[trajectoryProgress];
        } else {
            componentCounts[reactant]--;
        }
    }

    for (i = 0; i < reactions[j].products.size(); i++) {
        if (!componentHasTrajectory(reactions[j].products[i].index)) {
            componentCounts[reactions[j].products[i].index] += reactions[j].products[i].change;
        }
    }
}

void Simulation::printReactions() {
    int i, j;

    printf("\nThe following reactions are simulated:\n\n");
    for (i = 0; i < sys.numReactions; i++) {
        if (reactions[i].reactants.empty())
            printf("0");
        else
            printf("%s ", (*componentNames)[reactions[i].reactants[0].index].c_str());
        for (j = 1; j < reactions[i].reactants.size(); j++)
            printf("+ %s ", (*componentNames)[reactions[i].reactants[j].index].c_str());
        printf(" ->  ");
        if (reactions[i].products.empty())
            printf("0 ");
        else
            printf("%2d %s ", reactions[i].products[0].change,
                   (*componentNames)[reactions[i].products[0].index].c_str());
        for (j = 1; j < reactions[i].products.size(); j++)
            printf("+ %2d %s ", reactions[i].products[j].change,
                   (*componentNames)[reactions[i].products[j].index].c_str());
        printf("\t\tk = %4.5f\n", reactions[i].reactionConstant);
    }
}

void Simulation::printTrajectory(std::ostream &os, Trajectory &trajectory) {
    auto jsonObj = trajectory.getJson();
    auto reactionsJson = nlohmann::json();
    for (auto r : reactions) {
        auto rJson = nlohmann::json();
        auto reactants = std::vector<std::string>();
        for (auto react : r.reactants) {
            reactants.push_back((*componentNames)[react.index]);
        }
        rJson["reactants"] = nlohmann::json(reactants);
        rJson["k"] = nlohmann::json(r.reactionConstant);
        reactionsJson.push_back(rJson);
    }
    jsonObj["reactions"] = reactionsJson;
    jsonObj["random_variates"] = randomVariates;

    nlohmann::json::to_msgpack(jsonObj, os);
}

/**
 *
 * @return The number of components used in this simulation.
 */
unsigned long Simulation::getNumComponents() {
    return sys.numComponents;
}
