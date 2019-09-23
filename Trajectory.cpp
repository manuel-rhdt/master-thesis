//
// Created by reinhardt on 17-09-19.
//

#include <cassert>
#include <iomanip>
#include <limits>
#include <iostream>
#include "Trajectory.hh"

Trajectory::Trajectory(unsigned long numComponents) : timeStamps(vector<double>()) {
    componentCounts = vector<vector<double>>(numComponents, vector<double>());
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

    for (auto &componentArray : componentCounts) {
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

    return componentCounts[component][timeIndex];
}

std::ostream &operator<<(std::ostream &os, const Trajectory &trajectory) {
    // Save output stream state, so we can reset if at the end of this function
    std::ios_base::fmtflags flags(os.flags());
    os << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
    os << "# Trajectory" << endl;
    os << "# " << trajectory.timeStamps.size() << " timestamps; " << trajectory.componentCounts.size()
       << " components." << endl;
    for (auto t : trajectory.timeStamps) {
        os << t << " ";
    }
    os << endl;
    // for the component concentrations we just want to print integers (TODO: figure out whether to uncomment the following line)
//    os << std::fixed << std::setprecision(0);
    for (const auto &componentArray : trajectory.componentCounts) {
        for (auto c : componentArray) {
            os << c << " ";
        }
        os << endl;
    }

    os.flags(flags);

    return os;
}


std::istream &operator>>(std::istream &is, Trajectory &trajectory) {
    vector<vector<double>> numArray = vector<vector<double >>();
    numArray.emplace_back(vector<double>());
    is >> std::noskipws;
    while (!is.eof()) {
        double num = 0.0;
        if ((is >> num).fail()) {
            // we either reach end of line or '#' which we interpret as ignoring the rest of the line.
            is.clear();

            char failedChar;
            is >> failedChar;
            if (failedChar == ' ' || failedChar == '\t') {
                continue;
            } else if (failedChar == '\n' || failedChar == '#') {
                if (failedChar != '\n') {
                    is.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                if (!numArray.back().empty()) {
                    numArray.emplace_back(vector<double>());
                }
                continue;
            } else {
                if (is.eof()) {
                    break;
                }
                std::cerr << "Found invalid character '" << (int) failedChar << "' while parsing trajectory.";
                std::terminate();
            }
        }
        if (is) {
            numArray.back().push_back(num);
        }
    }
    if (numArray.back().empty()) {
        numArray.pop_back();
    }

    trajectory.timeStamps = numArray.front();
    numArray.erase(numArray.begin());
    trajectory.componentCounts = std::move(numArray);

    is.clear();
    return is;
}

/**
 * Construct an empty trajectory.
 *
 * This is useful if you want to read in a trajectory from a file.
 */
Trajectory::Trajectory() {

}


double TrajectoryIterator::next() {

}