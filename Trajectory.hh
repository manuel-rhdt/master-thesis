//
// Created by reinhardt on 17-09-19.
//

#ifndef GILLESPIE_TRAJECTORY_HH
#define GILLESPIE_TRAJECTORY_HH

#include <vector>
#include <ostream>
#include <json/single_include/nlohmann/json.hpp>

using namespace std;


class Trajectory {
public:
    Trajectory();

    explicit Trajectory(unsigned long numComponents);

    /**
     * A vector of vectors to store the concentration changes for every component.
     *
     * `componentConcentrations` is a list that contains `numComponents` elements. Every element itself is a list with
     * the same length as `timeStamps` and contains the corresponding component concentration at each time stamp.
     */
    vector<vector<double>> componentCounts{vector<vector<double>>()};

    vector<double> timeStamps{vector<double>()};

    vector<unsigned int> reactions{vector<unsigned int>()};

    void insertTimestamp(double dt);

    // read/write trajectories from/to streams in an ASCII format (e.g. files)

    friend std::ostream &operator<<(std::ostream &os, const Trajectory &trajectory);

    friend std::istream &operator>>(std::istream &is, Trajectory &trajectory);

    [[nodiscard]] nlohmann::json getJson() const;
};


#endif //GILLESPIE_TRAJECTORY_HH
