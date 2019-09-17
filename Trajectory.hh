//
// Created by reinhardt on 17-09-19.
//

#ifndef GILLESPIE_TRAJECTORY_HH
#define GILLESPIE_TRAJECTORY_HH

#include <vector>
#include <ostream>

using namespace std;


class Trajectory {
private:

    vector<double> timeStamps{vector<double>()};
public:
    explicit Trajectory(unsigned long numComponents);

    /**
     * A vector of vectors to store the concentration changes for every component.
     *
     * `componentConcentrations` is a list that contains `numComponents` elements. Every element itself is a list with
     * the same length as `timeStamps` and contains the corresponding component concentration at each time stamp.
     */
    vector<vector<double>> componentConcentrations{vector<vector<double>>()};

    void insertTimestamp(double dt);

    double getComponentConcentrationAt(double time, unsigned long component);

    friend std::ostream &operator<<(std::ostream &os, const Trajectory &trajectory);
};


#endif //GILLESPIE_TRAJECTORY_HH
