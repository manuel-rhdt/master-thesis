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
    int runType;
    double elapsedTime;
    std::vector<statistics> stat;
    const std::vector<std::string> *componentNames;

public:
    Run(int runType, std::vector<std::string> *componentNames);

    ~Run();
};

#endif //GILLESPIE_RUN_HH
