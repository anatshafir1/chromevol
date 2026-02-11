#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <unordered_map>

using namespace bpp;
using namespace std;

// 1. Need to test likelihood calculation before the actual implementation
// 2. Need to think whether this implementation will handle other models, such as OU, for example (given a simple transformation)

// The algorithm:
// 1. create an inordered map of each node id and the respective branch. This map will be
// updated throughout the algorithm
double calculateSigma(vector<double> &u, vector<double> &V){
    double sigma = 0;
    size_t n = u.size();
    for (size_t i = 0; i < n; i++){
        sigma += ((u[i]*u[i])/(V[i]*V[i]));
    }
    return sigma/static_cast<double>(n);
}
double getExpectedTraitValue(double &t1, double &t2, double &y1, double &y2){
    double nominator = ((1/t1)*y1) + ((1/t2)*y2);
    double denominator = (1/t1) + (1/t2);
    return nominator/denominator;

}
double calculateNewBranchLength(double& originalBranchLength, double &t1, double &t2){
    double newBranchLength = originalBranchLength;
    newBranchLength += ((t1*t2)/(t1+t2));
    return newBranchLength;
}

void createNodeIdsTraitValuesMap(shared_ptr<PhyloTree> tree, std::unordered_map<uint, double> &traitValues, std::unordered_map<string, double> &leafTraitStates){
    auto leaves = tree->getAllLeaves();
    string leafName;
    for (auto &leaf : leaves){
        leafName = leaf->getName();
        traitValues[tree->getNodeIndex(leaf)] = leafTraitStates[leafName];
    }
}
void createBranchLengthsMap(shared_ptr<PhyloTree> tree, std::unordered_map<uint, double> &branchLengths){
    //auto father = tree_->getFatherOfNode(tree_->getNode(nodeId));
    //uint fatherId = tree_->getNodeIndex(father);
    auto nodeIds = tree->getNodeIndexes(tree->getAllNodes());
    for (auto &nodeId : nodeIds){
        if (nodeId == tree->getRootIndex()){
            continue;
        }
        if (tree->isLeaf(nodeId)){
            auto branch = tree->getEdgeToFather(nodeId);
            double branchLength = branch->getLength();
            branchLengths[nodeId] = branchLength;

        }else{
            branchLengths[nodeId] = -1; // should be updated during the recursion
        }

    }
}

void updateRequiredVectorsForLikCalculationRec(shared_ptr<PhyloTree> tree, uint nodeId, std::unordered_map<uint, double> &branchLengths, std::unordered_map<uint, double> &traitValues, vector<double> &u, vector<double> &V){
    if (tree->isLeaf(nodeId)){
        return;
    }
    auto sons = tree->getSons(nodeId);
    vector<double> sonsBranchLengths;
    for (auto &sonId : sons){
        if (!(tree->isLeaf(sonId))){
            updateRequiredVectorsForLikCalculationRec(tree, sonId, branchLengths, traitValues, u, V);
           
        }
    }
    double t1 = branchLengths[sons[0]];
    double t2 = branchLengths[sons[1]];
    if ((t1 < 0) || (t2 < 0)){
        throw Exception("updateRequiredVectorsForLikCalculationRec(): Something went wrong in the recursion!");

    }
    u.push_back(t2-t1);
    V.push_back(t1+t2);
    traitValues[nodeId] = getExpectedTraitValue(t1, t2, traitValues[sons[0]], traitValues[sons[1]]);
    auto branch = tree->getEdgeToFather(nodeId);
    double branchLength = branch->getLength();
    branchLengths[nodeId] = branchLength + calculateNewBranchLength(branchLength, t1, t2);
    if (nodeId == tree->getRootIndex()){
        u.push_back(0);
        V.push_back(branchLengths[nodeId]);
    }
    
}

double calculateLikelihood(shared_ptr<PhyloTree> tree, std::unordered_map<string, double> &leafTraitStates){
    std::unordered_map<uint, double> branchLengths;
    std::unordered_map<uint, double> traitValues;
    createNodeIdsTraitValuesMap(tree, traitValues, leafTraitStates);
    createBranchLengthsMap(tree, branchLengths);
    vector<double> V; // variances (same notation as in the paper)
    vector<double> u; // contrasts (same notation as in the paper)
    uint rootId = tree->getRootIndex();
    updateRequiredVectorsForLikCalculationRec(tree, rootId, branchLengths, traitValues, u, V);
    double sigma = calculateSigma(u, V);
    double likelihood = 0;
    for (size_t i = 0; i < u.size(); i++){
        likelihood += (std::log(V[i]) + ((u[i]*u[i])/(sigma * V[i])));
    }
    likelihood += static_cast<double>(u.size())*std::log(2* M_PI);
    likelihood *= (-0.5);
    return likelihood;

}



int main(){
    Newick reader;
    shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("((((P:0.21,H:0.21):0.28,M:0.49):0.13,A:0.62):0.38,G:1.0);", false, "", false, false));
    std::unordered_map<string, double> traitData;
    traitData["G"] = -1.46968;
    traitData["A"] = 2.02815;
    traitData["M"] = 2.37024;
    traitData["P"] = 3.61092;
    traitData["H"] = 4.09434;
    double likelihood = calculateLikelihood(pTree, traitData);
    std::cout << "likelihood is " << likelihood << std::endl;


}