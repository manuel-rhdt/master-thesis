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

nlohmann::json Trajectory::getJson() const {
    auto output = nlohmann::json();
    output["timestamps"] = timeStamps;
    output["reaction_events"] = reactions;
    output["components"] = componentCounts;
    return output;
}

std::ostream &operator<<(std::ostream &os, const Trajectory &trajectory) {
    os << trajectory.getJson();
    return os;
}

std::istream &operator>>(std::istream &is, Trajectory &trajectory) {
    auto jsonObj = nlohmann::json();
    is >> jsonObj;

    trajectory.timeStamps = jsonObj["timestamps"].get<vector<double>>();
    trajectory.componentCounts = jsonObj["components"].get<vector<vector<double>>>();

    return is;
}

/**
 * Construct an empty trajectory.
 *
 * This is useful if you want to read in a trajectory from a file.
 */
Trajectory::Trajectory() = default;


