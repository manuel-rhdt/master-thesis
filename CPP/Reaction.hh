//
// Created by reinhardt on 20-09-19.
//

#ifndef GILLESPIE_REACTION_HH
#define GILLESPIE_REACTION_HH


#include <vector>

struct ReactionComponent {
    int index;
    int change;
};

class Reaction {
public:
    double reactionConstant{0.0};
    std::vector<ReactionComponent> reactants{std::vector<ReactionComponent>()};
    std::vector<ReactionComponent> products{std::vector<ReactionComponent>()};
};


#endif //GILLESPIE_REACTION_HH
