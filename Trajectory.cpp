//
// Created by reinhardt on 17-09-19.
//

#include <cassert>
#include <iomanip>
#include <limits>
#include <iostream>
#include "Trajectory.hh"
#include "nlohmann/json.hpp"

Trajectory::Trajectory(unsigned long numComponents) : timeStamps(vector<double>()) {
    componentCounts = vector<vector<double>>(numComponents, vector<double>());
}

/**
 * Inserts a new timestamp exactly `dt` after the last timestamp.
 *
 * The component concentrations will initially be zero for the new timestamp.
 */
void Trajectory::insertTimestamp(double dt) {
    double previousTimestamp = 0.0;
    if (timeStamps.size() > 0) {
        previousTimestamp = timeStamps.back();
    }
    timeStamps.push_back(previousTimestamp + dt);

    reactions.push_back(0);

    for (auto &componentArray : componentCounts) {
        componentArray.push_back(0.0);
    }
}

std::vector<double> &Trajectory::getComponent(std::string &name) {
    for (long comp = 0; comp < componentCounts.size(); comp++) {
        if ((*componentNames)[comp] == name) {
            return componentCounts[comp];
        }
    }
    throw std::runtime_error("No component with name " + name);
}

nlohmann::json Trajectory::getJson() const {
    if (!componentNames || componentNames->empty()) {
        throw std::runtime_error("The trajectory does not know the names of its components");
    }
    auto output = nlohmann::json();
    output["timestamps"] = timeStamps;
    output["reaction_events"] = reactions;
    for (int compNum = 0; compNum < componentCounts.size(); compNum++) {
        auto name = (*componentNames)[compNum];
        output["components"][name] = nlohmann::json(componentCounts[compNum]);
    }
    return output;
}

std::ostream &operator<<(std::ostream &os, const Trajectory &trajectory) {
    os << trajectory.getJson();
    return os;
}

std::istream &operator>>(std::istream &is, Trajectory &trajectory) {
    auto jsonObj = nlohmann::json::from_msgpack(is);

    trajectory.timeStamps = jsonObj["timestamps"].get<vector<double>>();

    trajectory.componentNames = make_shared<vector<string>>(vector<string>());
    trajectory.componentCounts = vector<vector<double>>();
    for (auto const &comp : jsonObj["components"].get<map<string, vector<double>>>()) {
        trajectory.componentNames->push_back(comp.first);
        trajectory.componentCounts.push_back(comp.second);
    }

    return is;
}

/**
 * Construct an empty trajectory.
 *
 * This is useful if you want to read in a trajectory from a file.
 */
Trajectory::Trajectory() = default;


