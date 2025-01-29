#include "BrownianMotionLikelihood.h"

using namespace std;
using namespace bpp;

void BrownianMotionLikelihood::calculateContrasts(){
    calculatedContrasts_ = true;
    size_t index = 0;
    updateContrastsVectors(index, tree_->getRootIndex());
    calculateSigmaMLE();
    calculateMuMLE();
}

double BrownianMotionLikelihood::getExpectedTraitValue(double &t1, double &t2, double &trait1, double &trait2){
    double nominator = ((1 / t1) * trait1) + ((1 / t2) * trait2);
    double denominator = (1 / t1) + (1 / t2);
    return nominator / denominator;

}
double BrownianMotionLikelihood::calculateNewBranchLength(double branchLength, double &t1, double &t2){
    return branchLength + ((t1 * t2) / (t1 + t2));

}
void BrownianMotionLikelihood::updateContrastsVectors(size_t &index, uint nodeId){
    if (tree_->isLeaf(tree_->getNode(nodeId))){
        auto branch = tree_->getEdgeToFather(nodeId);
        double branchLength = branch->getLength();
        correctedBranchLengths_[nodeId] = branchLength;
        return;
    }
    auto sons = tree_->getSons(nodeId);
    for (auto &sonId : sons){
        updateContrastsVectors(index, sonId);
    }
    u_[index] = parsimonyStates_[sons[0]] - parsimonyStates_[sons[1]];
    V_[index] = correctedBranchLengths_[sons[0]] + correctedBranchLengths_[sons[1]];
    index ++;
    double branchLength;
    if (nodeId == tree_->getRootIndex()){
        branchLength = 0;
    }else{
        auto branch = tree_->getEdgeToFather(nodeId);
        branchLength = branch->getLength();
    }
    parsimonyStates_[nodeId] = getExpectedTraitValue(correctedBranchLengths_[sons[0]], correctedBranchLengths_[sons[1]], parsimonyStates_[sons[0]], parsimonyStates_[sons[1]]);
    correctedBranchLengths_[nodeId] = calculateNewBranchLength(branchLength, correctedBranchLengths_[sons[0]], correctedBranchLengths_[sons[1]]);
    if (nodeId == tree_->getRootIndex()){
        u_[index] = 0;
        V_[index] = correctedBranchLengths_[nodeId];
    }
    
}