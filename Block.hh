//
//  Block.hh
//  Gillespie
//
//  Created by Manuel Reinhardt on 2019-09-16.
//  Copyright © 2019 Manuel Reinhardt. All rights reserved.
//

#ifndef GILLESPIE_BLOCK_HH
#define GILLESPIE_BLOCK_HH

#include "Gillespie.hh"
#include "Run.hh"
#include <vector>

class Block {
private:
    double elapsedTime;
    std::vector<statistics> stat;
    const std::shared_ptr<std::vector<std::string>> componentNames;
public:

    Block(unsigned long numComponents, std::shared_ptr<std::vector<std::string>> componentNames);

    void updateRun(Run &run);

    friend std::ostream &operator<<(std::ostream &os, const Block &obj);

    void accumulate(double dt, const std::vector<double> &componentCounts);
};

#endif /*GILLESPIE_BLOCK_HH*/
