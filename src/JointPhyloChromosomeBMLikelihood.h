//
// File: FormulaOfPhyloLikelihood.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: jeudi 8 dÃÂ©cembre 2016, ÃÂ  10h 35
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

#ifndef CHROMEVOL_JOINTPHYLOCHROMOSOMEBMLIKELIHOOD_H
#define CHROMEVOL_JOINTPHYLOCHROMOSOMEBMLIKELIHOOD_H

#include <Bpp/Numeric/Function/Operators/ComputationTree.h>

#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowNumeric.h>

#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SetOfAbstractPhyloLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include "BrownianMotionAncestralReconstruction.h"
#include "BrownianMotionLikelihood.h"
#include "ChromosomeSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The JointPhyloLikelihood class is intended for likelihood computation of the following form- given two types of data in the leaves,
 * D1 and D2, if the two evolutionary process are not independent, then if for example, D1= {0,1}, each state dictates a different evolutionary process for D2, therefore the formula is:
 * P(D1,D2| teta) = P(D1|teta)*P(D2|D1,teta). The first term can be calculated via the Felsenstein algorithm, while in order to compute the second term, we have to know the history
 * of the transitions of D1 anlong the phylogeny. This history is unknown, but can be estimated using a large number of stochastic mappings, that can be summarized to generate
 * an expected history of trajectories. Knowing the expected history will allow to assign each model of D2 according to the states of D1, and thus enable to calculate the second term using Felsenstein
 * 
 *
 * WARNING: This formula applies on the log-likelihoods (ie getValues())
 *
 */ 

class JointPhyloChromosomeBMLikelihood :
  public SetOfAbstractPhyloLikelihood

{
protected:

  std::shared_ptr<LikelihoodCalculation> likCal_;
  bool muOptimization_;
  bool traitParamsOptimization_;
  bool updateAncestralsInModels_;
  std::shared_ptr<PhyloTree> tree_;
  std::shared_ptr<BrownianMotionLikelihood> bmLikelihood_;
  const std::unordered_map<string, double> continuousTraitStates_;
  std::shared_ptr<NumericMutable<DataLik>> bmLikNode_;
  bool setBMLikNode_;
  double minTraitState_;
  double maxTraitState_;


public:
  JointPhyloChromosomeBMLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, std::shared_ptr<LikelihoodCalculationSingleProcess> likCalculationProcess, std::shared_ptr<PhyloTree> tree, const std::unordered_map<string, double> &traitData, bool inCollection = true):
  AbstractPhyloLikelihood(context),
  SetOfAbstractPhyloLikelihood(context, pC, {}, inCollection),
  likCal_(new LikelihoodCalculation(context)),
  muOptimization_(false),
  traitParamsOptimization_(false),
  updateAncestralsInModels_(false),
  tree_(tree),
  bmLikelihood_(),
  continuousTraitStates_(traitData),
  bmLikNode_(),
  setBMLikNode_(false),
  minTraitState_(0),
  maxTraitState_(0)
  {
    addPhyloLikelihood(1, "");
    bmLikelihood_ = std::make_shared<BrownianMotionLikelihood>(tree, traitData);
    likCal_->setLikelihoodNode(makeLikelihoods());
    auto bmModel = std::dynamic_pointer_cast<const ChromosomeBMSubstitutionModel>((dynamic_cast<SingleProcessPhyloLikelihood*>(getAbstractPhyloLikelihood(nPhylo_[0]))->getSubstitutionProcess()).getModel(1));
    minTraitState_ = bmModel->getMinTraitState();
    maxTraitState_ = bmModel->getMaxTraitState();
    

  }

  ~JointPhyloChromosomeBMLikelihood() {
    SingleProcessPhyloLikelihood* lik = dynamic_cast<SingleProcessPhyloLikelihood*>(getAbstractPhyloLikelihood(nPhylo_[0]));
    auto sequenceDataChromosome = lik->getData();
    auto processChromosome= &(lik->getSubstitutionProcess());
    auto context = &(lik->getContext());
    delete processChromosome;
    delete sequenceDataChromosome;
    delete context;

  }

  JointPhyloChromosomeBMLikelihood* clone() const
  {
    return new JointPhyloChromosomeBMLikelihood(*this);
  }


  JointPhyloChromosomeBMLikelihood(const JointPhyloChromosomeBMLikelihood& sd) : 
    AbstractPhyloLikelihood(sd),
    SetOfAbstractPhyloLikelihood(sd),
    likCal_(sd.likCal_),
    muOptimization_(sd.muOptimization_),
    traitParamsOptimization_(sd.traitParamsOptimization_),
    updateAncestralsInModels_(sd.updateAncestralsInModels_),
    tree_(sd.tree_),
    bmLikelihood_(sd.bmLikelihood_),
    continuousTraitStates_(sd.continuousTraitStates_),
    setBMLikNode_(sd.setBMLikNode_),
    minTraitState_(sd.minTraitState_),
    maxTraitState_(sd.maxTraitState_)

{
  setLikelihoodBMNode();

}

public:


  /**
   * @brief Get the logarithm of the likelihood for the whole dataset.
   *
   * @return The logarithm of the likelihood of the dataset.
   */
  std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation() const
  {
    return likCal_;
  }

  void fireParameterChanged(const ParameterList& params);
  void setMuOptimization(){
    muOptimization_ = true;
  }
  void unsetMuOptimization(){
    muOptimization_ = false;
  }
  void setTraitOptimization(){
    traitParamsOptimization_ = true;
  }
  void setChromosomeOptimization(){
    traitParamsOptimization_ = false;
  }
  const double getMinTraitState() const{return minTraitState_;}
  const double getMaxTraitState() const{return maxTraitState_;}
  void getInitialLikelihood();
  const std::shared_ptr<BrownianMotionLikelihood> getBMLikelihoodObject() const {return bmLikelihood_;}

  



protected:
  /**
   * @brief Build the LikelihoodNode from the computation Tree
   *
   */
  ValueRef<DataLik> makeLikelihoods();
  void setLikelihoodBMNode();
  void updateAllAncestralStatesInBranches();




};
} // end of namespace bpp.
#endif // CHROMEVOL_JOINTPHYLOCHROMOSOMEBMLIKELIHOOD_H
