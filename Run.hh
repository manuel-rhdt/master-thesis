//
// Created by reinhardt on 13-09-19.
//

#ifndef GILLESPIE_RUN_HH
#define GILLESPIE_RUN_HH

#include "Gillespie.hh"
#include <string>
#include <vector>

class Run {
private:
    const std::vector<std::string> componentNames;

public:
    Run(int runType, std::vector<std::string> componentNames);

    Run(const Run &) = delete;

    void finish(void);

    std::vector<statistics> stat;
    double elapsedTime;
};

#endif //GILLESPIE_RUN_HH
