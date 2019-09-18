//
// Created by reinhardt on 17-09-19.
//

#include <cassert>
#include <iomanip>
#include <limits>
#include "Trajectory.hh"

Trajectory::Trajectory(unsigned long numComponents) : timeStamps(vector<double>()) {
    componentConcentrations = vector<vector<double>>(numComponents, vector<double>());
}

/**
 * Inserts a new timestamp exactly `dt` after the last timestamp.
 */
void Trajectory::insertTimestamp(double dt) {
    double previousTimestamp = 0.0;
    if (timeStamps.size() > 0) {
        previousTimestamp = timeStamps.back();
    }
    timeStamps.push_back(previousTimestamp + dt);

    for (auto &componentArray : componentConcentrations) {
        componentArray.push_back(0.0);
    }
}

double Trajectory::getComponentConcentrationAt(double time, unsigned long component) {
    int timeIndex = 0;
    while (timeStamps[timeIndex] < time) {
        timeIndex++;
    }
    // we want the concentrations just before `time`
    timeIndex--;

    assert(timeIndex > 0);

    return componentConcentrations[component][timeIndex];
}

std::ostream &operator<<(std::ostream &os, const Trajectory &trajectory) {
    // Save output stream state, so we can reset if at the end of this function
    std::ios_base::fmtflags flags(os.flags());
    os << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
    os << "# Trajectory" << endl;
    os << "# " << trajectory.timeStamps.size() << " timestamps; " << trajectory.componentConcentrations.size()
       << " components." << endl;
    for (auto t : trajectory.timeStamps) {
        os << t << " ";
    }
    os << endl;
    // for the component concentrations we just want to print integers
    os << std::fixed << std::setprecision(0);
    for (const auto &componentArray : trajectory.componentConcentrations) {
        for (auto c : componentArray) {
            os << c << " ";
        }
        os << endl;
    }

    os.flags(flags);

    return os;
}