//
//  Block.cpp
//  Gillespie
//
//  Created by Manuel Reinhardt on 2019-09-16.
//  Copyright © 2019 Manuel Reinhardt. All rights reserved.
//

#include <iomanip>
#include <memory>
#include "Block.hh"
#include "Run.hh"

statistics getAverageStatistics(statistics item, double elapsedTime) {
    double sum = item.sum / elapsedTime;
    double sumSquared = item.sumsq / elapsedTime;
    double sigmaSquared = sumSquared - sum * sum; /*sigma^2*/
    double err = 0.0;
    double noise = 0.0;
    if (sigmaSquared > 0.0) {
        err = sqrt(sigmaSquared); /*sigma*/
    }

    if (sum != 0.0) noise = err / sum;

    return statistics{
            item.acc, sum, sumSquared, err, noise
    };
}

Block::Block(std::vector<statistics>::size_type numComponents,
             const std::shared_ptr<std::vector<std::string>> componentNames)
        : elapsedTime(0.0), componentNames(componentNames) {
    stat = std::vector<statistics>(numComponents, statistics{0, 0.0, 0.0, 0.0, 0.0});
}

void Block::accumulate(double dt, const std::vector<double> &componentCounts) {
    elapsedTime += dt;
    for (unsigned long i = 0; i < stat.size(); i++) {
        stat[i].sum += componentCounts[i] * dt;
        stat[i].sumsq += componentCounts[i] * componentCounts[i] * dt;
        stat[i].acc++;
    }
}

void Block::updateRun(Run &run) {
    run.elapsedTime += elapsedTime;
    for (std::vector<statistics>::size_type i = 0; i < stat.size(); i++) {
        statistics runStat = getAverageStatistics(stat[i], elapsedTime);
        run.stat[i].sum += runStat.sum;
        run.stat[i].sumsq += runStat.sum * runStat.sum;
        run.stat[i].noise += runStat.noise;
        run.stat[i].acc++;
    }
}

/**
 * Write out the block
 *
 * @param os
 * @param obj
 * @return
 */
std::ostream &operator<<(std::ostream &os, const Block &obj) {
    os << "Block\n";
    os << std::string("Component\tblock average\t   deviation\t       noise\t     #counts\n");
    for (unsigned long comp = 0; comp < obj.stat.size(); comp++) {
        statistics runStat = getAverageStatistics(obj.stat[comp], obj.elapsedTime);
        os << std::scientific << std::setprecision(4) << std::setw(9)
           << (*obj.componentNames)[comp] << "\t" << std::setw(13)
           << runStat.sum / obj.elapsedTime << '\t' << std::setw(12)
           << runStat.err << '\t' << std::setw(12)
           << runStat.noise << '\t' << std::setw(12)
           << runStat.acc << '\n';
    }

    return os;
}