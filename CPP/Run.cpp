//
// Created by reinhardt on 13-09-19.
//

#include "Run.hh"
#include <cassert>

Run::Run(int numComponents, const std::vector<std::string> componentNames)
        : componentNames(componentNames), elapsedTime(0.0) {
    stat = std::vector<statistics>(numComponents, statistics{0, 0.0, 0.0, 0.0, 0.0});
}

void Run::finish(void) {
    printf("\n==============================================================================\n\n");

    printf("Elapsed time         %8.4f\n\n", elapsedTime);

    printf("Run averages.\n\n");
    printf("Component\tRun average\t   error\t   noise\t\t#counts\n");

    for (auto item = stat.begin(); item != stat.end(); ++item) {
        item->sum /= item->acc;
        item->sumsq /= item->acc;
        item->err = (item->sumsq - item->sum * item->sum) / (double) item->acc;
        if (item->err > 0.)
            item->err = sqrt(item->err);
        item->noise /= item->acc;
        printf("%s\t\t%8.4f\t%8.4f\t%8.4f\t%8d\n", componentNames[item - stat.begin()].c_str(), item->sum, item->err,
               item->noise, item->acc);
    }

    stat.clear();
}