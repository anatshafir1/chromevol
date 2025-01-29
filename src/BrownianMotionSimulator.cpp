#include "BrownianMotionSimulator.h"

using namespace std;
using namespace bpp;

void BrownianMotionSimulator::simulate(){
    std::default_random_engine generator;
    uint rootId = tree_->getRootIndex();
    minStateAbsolute_ = mu_;
    maxStateAbsolute_ = mu_;
    maxStateInLeaves_ = -std::numeric_limits<double>::infinity();
    minStateInLeaves_ = std::numeric_limits<double>::infinity();
    traitStatesInNodes_[rootId] = mu_;
    generator.seed(seed_);
    simulate_(rootId, generator);
}
void BrownianMotionSimulator::simulate_(uint &nodeId, std::default_random_engine &generator){
    if (tree_->isLeaf(nodeId)){
        if (traitStatesInNodes_.at(nodeId) > maxStateInLeaves_){
            maxStateInLeaves_ = traitStatesInNodes_.at(nodeId);
        }else if (traitStatesInNodes_.at(nodeId) < minStateInLeaves_){
            minStateInLeaves_ = traitStatesInNodes_.at(nodeId);
        }
        return;
    }
    auto sons = tree_->getSons(nodeId);
    for (auto &sonId: sons){
        auto branch = tree_->getEdgeToFather(sonId);
        double branchLength = branch->getLength();
        double standard_deviation = std::sqrt(sigma_ * branchLength);
        std::normal_distribution<double> distribution(0.0, standard_deviation);
        double x = distribution(generator);
        double sonState = x + traitStatesInNodes_.at(nodeId);
        if (sonState > maxStateAbsolute_){
            maxStateAbsolute_ = sonState;
        }else if (sonState < minStateAbsolute_){
            minStateAbsolute_ = sonState;
        }
        traitStatesInNodes_[sonId] = sonState;
        traitStatesBranches_[sonId] = (traitStatesInNodes_.at(nodeId) + sonState)/2;
        simulate_(sonId, generator);
    }
}

std::unordered_map<string, double> BrownianMotionSimulator::getLeavesStates() const{
    std::unordered_map<string, double> leavesStates;
    auto it = traitStatesInNodes_.begin();
    while (it != traitStatesInNodes_.end()){
        auto nodeId = it->first;
        if (tree_->isLeaf(nodeId)){
            leavesStates[tree_->getNode(nodeId)->getName()] = traitStatesInNodes_.at(nodeId);
        }
        it ++;
    }
    return leavesStates;

}