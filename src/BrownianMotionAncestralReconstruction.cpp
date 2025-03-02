//
// File: BrownianMotionAncestralReconstruction.cpp
// Authors:
//   Anat Shafir
// Created: 2024-12-11 15:09:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include "BrownianMotionAncestralReconstruction.h"

using namespace std;
using namespace bpp;

void BrownianMotionAncestralReconstruction::reconstructAncestralStates(){
    addLeavesToAncestralStates();
    uint rootId = tree_->getRootIndex();
    reconstructAncestralStatesPostOrder(rootId);
    reconstructAncestralStatesPreOrder(rootId);
    //revertLeafStates();
    correctMLStatesAccordingToMu();
    
}
void BrownianMotionAncestralReconstruction::claculateNodeDepthsRec(std::unordered_map<uint, double> &nodeDepths, uint nodeId, uint sonId){
    
    if (nodeDepths[nodeId] < 0){
        nodeDepths[nodeId] = nodeDepths[sonId];
        auto branch = tree_->getEdgeToFather(sonId);
        nodeDepths[nodeId] = nodeDepths[sonId] + branch->getLength();
    }else{
        return;
    }
    if (nodeId == tree_->getRootIndex()){
        return;
    }
    auto fatherNode = tree_->getFatherOfNode(tree_->getNode(nodeId));
    uint fatherIndex = tree_->getNodeIndex(fatherNode);
    claculateNodeDepthsRec(nodeDepths, fatherIndex, nodeId);

}
void BrownianMotionAncestralReconstruction::claculateNodeDepths(std::unordered_map<uint, double> &nodeDepths){
    auto nodes = tree_->getAllNodes();
    for (auto &node : nodes){
        nodeDepths[tree_->getNodeIndex(node)] = -1;
    }
    auto leafNodes = tree_->getLeavesUnderNode(tree_->getNode(tree_->getRootIndex()));
    for (auto &leaf : leafNodes){
        uint nodeId = tree_->getNodeIndex(leaf);
        nodeDepths[nodeId] = 0;
        auto fatherNode = tree_->getFatherOfNode(tree_->getNode(nodeId));
        uint fatherIndex = tree_->getNodeIndex(fatherNode);
        claculateNodeDepthsRec(nodeDepths, fatherIndex, nodeId);
    }

}
void BrownianMotionAncestralReconstruction::correctMLStatesAccordingToMu(){
    if (mu_ == muMLE_){
        return;
    }
    std::unordered_map<uint, double> nodeDepths;
    claculateNodeDepths(nodeDepths);
    double rootDepth = nodeDepths[tree_->getRootIndex()];
    auto nodeIds = tree_->getNodeIndexes(tree_->getAllNodes());
    for (auto &nodeId : nodeIds){
        ancestralTraitStates_[nodeId] += (mu_ - muMLE_) * (nodeDepths[nodeId] / rootDepth);

    }
}

// void BrownianMotionAncestralReconstruction::revertLeafStates(){
//     addLeavesToAncestralStates();

// }
void BrownianMotionAncestralReconstruction::addLeavesToAncestralStates(){
    auto leafNodes = tree_->getLeavesUnderNode(tree_->getNode(tree_->getRootIndex()));
    for (auto &leaf : leafNodes){
        auto leafName = leaf->getName();
        uint nodeId = tree_->getNodeIndex(leaf);
        ancestralTraitStates_[nodeId] = leafTraitStates_[leafName];
    }
}

void BrownianMotionAncestralReconstruction::reconstructAncestralStatesPostOrder(uint nodeId){
    if (tree_->isLeaf(tree_->getNode(nodeId))){
        auto branch = tree_->getEdgeToFather(nodeId);
        double branchLength = branch->getLength();
        p_[nodeId] = 1/branchLength;
        //ancestralTraitStates_[nodeId] -= (muMLE_ - mu_); 

    }else{
        double p_A = 0;
        auto sons = tree_->getSons(nodeId);
        for (size_t i = 0; i < sons.size(); i++){
            reconstructAncestralStatesPostOrder(sons[i]);
            ancestralTraitStates_[nodeId] += (p_[sons[i]] * ancestralTraitStates_[sons[i]]);
            p_A += p_[sons[i]];
        }
        double branchLength;
        if (nodeId == tree_->getRootIndex()){
            branchLength = 0;
        }else{
            auto branch = tree_->getEdgeToFather(nodeId);
            branchLength = branch->getLength();
            
        }
        p_[nodeId] = p_A / (1 + branchLength * p_A);
        ancestralTraitStates_[nodeId] /= p_A;
        
    }
}

void BrownianMotionAncestralReconstruction::reconstructAncestralStatesPreOrder(uint nodeId){
    if(tree_->getRootIndex() != nodeId){
        auto fatherNode = tree_->getFatherOfNode(tree_->getNode(nodeId));
        uint fatherIndex = tree_->getNodeIndex(fatherNode);
        auto branch = tree_->getEdgeToFather(nodeId);
        double branchLength = branch->getLength();
        if (tree_->isLeaf(tree_->getNode(nodeId))){
            p_[nodeId] = 0;
            return;
        }else{
            ancestralTraitStates_[nodeId] = (ancestralTraitStates_[nodeId] * p_[nodeId] * branchLength) +
             ancestralTraitStates_[fatherIndex] - (ancestralTraitStates_[fatherIndex] * p_[nodeId] * branchLength);
            p_[nodeId] = (p_[nodeId] / (1 - branchLength * p_[nodeId]) +
                    (p_[fatherIndex] - p_[nodeId]) / (1 + branchLength * (p_[fatherIndex] - p_[nodeId])));

        }       
    }
    auto sons = tree_->getSons(nodeId);
    for (size_t i = 0; i < sons.size(); i++){
        reconstructAncestralStatesPreOrder(sons[i]);

    }

}
