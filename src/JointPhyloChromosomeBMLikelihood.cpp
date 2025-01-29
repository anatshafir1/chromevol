//
// File: JointPhyloLikelihood.cpp
// Authors:
//   Anat Shafir
// Created: 2023-04-16 11:57:00
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

#include "JointPhyloChromosomeBMLikelihood.h"

using namespace std;
using namespace bpp;


void JointPhyloChromosomeBMLikelihood::fireParameterChanged(const ParameterList& params)
  {
    getAbstractPhyloLikelihood(nPhylo_[0])->matchParametersValues(params);
    if (traitParamsOptimization_){
      if (updateAncestralsInModels_){
        // no need to update the likelihood each time, as eventually after updating all the branch states, the likelihood will be updated
        return;
      }
      double sigma = (getAbstractPhyloLikelihood(nPhylo_[0])->getParameters()).getParameter("Chromosome.sigma_1").getValue();
      double mu = (getAbstractPhyloLikelihood(nPhylo_[0])->getParameters()).getParameter("Chromosome.mu_1").getValue();
      bmLikelihood_->calculateLikelihood(mu, sigma);
      if (muOptimization_){
        updateAllAncestralStatesInBranches();

      }
    }else{
      if (updateAncestralsInModels_){
        return;
      }
    }
    likCal_->setLikelihoodNode(makeLikelihoods());
    getAbstractPhyloLikelihood(nPhylo_[0])->getLikelihoodNode();
  }

  void JointPhyloChromosomeBMLikelihood::getInitialLikelihood(){
    updateAllAncestralStatesInBranches();
    likCal_->setLikelihoodNode(makeLikelihoods());
    getAbstractPhyloLikelihood(nPhylo_[0])->getLikelihoodNode();

  }
  
  void JointPhyloChromosomeBMLikelihood::updateAllAncestralStatesInBranches(){
    // mu optimization induces the generation of new ML ancestral states in BM.


    double sigma = (getAbstractPhyloLikelihood(nPhylo_[0])->getParameters()).getParameter("Chromosome.sigma_1").getValue();
    double mu = (getAbstractPhyloLikelihood(nPhylo_[0])->getParameters()).getParameter("Chromosome.mu_1").getValue();
    bmLikelihood_->calculateLikelihood(mu, sigma);
    double sigmaMLE = bmLikelihood_->getSigmaMLE();
    double muMLE = bmLikelihood_->getMuMLE();
    std::shared_ptr<BrownianMotionAncestralReconstruction> bmAncestral = std::make_shared<BrownianMotionAncestralReconstruction>(muMLE, sigmaMLE, mu, sigma, continuousTraitStates_, tree_);
    bmAncestral->reconstructAncestralStates();
    auto reconstructedStates = bmAncestral->getAncestralStates(); // a map of node ids and the respective trait states
    // changing the state along each branch (the average between the two terminals of the branch):
    // 1. Calculate the average between the terminal states of each branch
    // 2. Need to get the id assigned for each model
    vector<uint> nodeIds = tree_->getNodeIndexes(tree_->getAllNodes());
    updateAncestralsInModels_ = true;
    for (auto &nodeId : nodeIds){
      if (nodeId == tree_->getRootIndex()){
        continue;
      }
      auto fatherNode = tree_->getFatherOfNode(tree_->getNode(nodeId));
      uint fatherIndex = tree_->getNodeIndex(fatherNode);
      double fatherState = reconstructedStates.at(fatherIndex);
      double currentNodeState = reconstructedStates.at(nodeId);
      double branchState = (currentNodeState + fatherState)/2;
      uint modelNum = static_cast<uint>((dynamic_cast<SingleProcessPhyloLikelihood*>(getAbstractPhyloLikelihood(nPhylo_[0]))->getSubstitutionProcess()).getModelNumberForNode(nodeId));
      setParameterValue("Chromosome.state_" + std::to_string(modelNum), branchState);
    }
    updateAncestralsInModels_ = false;
  }

 ValueRef<DataLik> JointPhyloChromosomeBMLikelihood::makeLikelihoods(){
    auto lik = getAbstractPhyloLikelihood(nPhylo_[0])->getLikelihoodNode();
    setLikelihoodBMNode();
    // add the likelihoods because the target values are log likelihoods
    auto jointLik = CWiseAdd<DataLik, std::tuple<DataLik, DataLik> >::create(context_, {lik, bmLikNode_}, Dimension<DataLik> ());
    return jointLik;

 }
 void JointPhyloChromosomeBMLikelihood::setLikelihoodBMNode(){
  double bmLik = bmLikelihood_->getLikelihood();
  ExtendedFloat bmLikEf = ExtendedFloat{bmLik,0};
  bmLikEf.normalize();
  if (!setBMLikNode_){
    bmLikNode_ = NumericMutable<DataLik>::create(context_, bmLikEf);
    setBMLikNode_ = true;
  }else{
    bmLikNode_->setValue(bmLikEf);
    
  }
 }
