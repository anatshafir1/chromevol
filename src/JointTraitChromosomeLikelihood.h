//
// File: JointTraitChromosomeLikelihood.h
// Authors:
//   Anat Shafir
// Created: 17/04/2023
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

#ifndef CHROMEVOL_JOINTTRAITCHROMOSOMELIKELIHOOD_H
#define CHROMEVOL_JOINTTRAITCHROMOSOMELIKELIHOOD_H

#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/JointPhyloLikelihood.h>
//#include "BaseNumberOptimizer.h"
//#include "ChromosomeNumberOptimizer.h"
#include "ChromosomeSubstitutionModel.h"
#include "ChromosomeNumberOptimizer.h"
#include "LikelihoodUtils.h"
#include "TraitOptimizer.h"
namespace bpp
{
/**
 * @brief The JointPhyloLikelihood class, for phylogenetic
 * likelihood on several independent data.
 *
 * WARNING: This formula applies on the log-likelihoods (ie getValues())
 *
 */ 

class JointTraitChromosomeLikelihood :
  public JointPhyloLikelihood,
  public TraitOptimizer
{
  private:
  std::string traitModel_; // Binary for binary models. Integer for models with more then two states (not implemented yet)
  std::map<uint, vector<int>> fixedParamsChromosome_;
  std::vector<string> fixedParamsTrait_;
  std::vector<int> chromosomeRateChangeType_;
  std::vector <unsigned int> baseNumberCandidates_;
  VectorSiteContainer* chromosomeVsc_;
  std::map<int, std::vector<std::pair<uint, int>>> sharedParamsChromosome_;
  


            // vector <SingleProcessPhyloLikelihood*> vectorOfLikelohoods_;
            // //vector <Context> vectorOfContexts_;
            // const PhyloTree* tree_;
            // const ChromosomeAlphabet* alphabet_;
            // const VectorSiteContainer* vsc_;
            // bool optimizeBaseNumber_;
            // vector<unsigned int> numOfPoints_;
            // vector<unsigned int> numOfIterations_;
            // vector<unsigned int> numOfPointsNextRounds_;
            // vector<unsigned int> numOfIterationsNextRounds_;
            // string typeOfOptimizer_;
            // string baseNumOptimizationMethod_;
            // mutable std::map<uint, uint> baseNumberUpperBound_;
            // double tolerance_;
            // bool standardOptimization_;
            // int BrentBracketing_;
            // vector <double> probsForMixedOptimization_;
            // std::map<uint, vector<int>> fixedParams_;
            // mutable std::map<int, std::vector<std::pair<uint, int>>> sharedParams_;
            // uint numOfShiftsForward_;
            // bool backwardPhaseStarted_;
            // // number of shifts. For each shift, a new added node
            // std::map<uint, vector<uint>> prevModelsPartitions_;
            // std::map<uint, std::vector<std::pair<string, double>>> prevModelParams_;
            // std::map<uint, std::pair<double, double>> prevModelsAICcLikValues_;
            // std::map<uint, std::vector<double>> prevModelsRootFrequencies_;



public:
  JointTraitChromosomeLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool expectedHistory, bool weightedFrequencies, size_t numOfMappings, std::string traitModel, std::vector<int> &rateChangeType, VectorSiteContainer* chromosomeVsc, std::vector<unsigned int> &baseNumberCandidates, bool inCollection = true);

  ~JointTraitChromosomeLikelihood() {
    //delete chromosomeVsc_;
  }

  JointTraitChromosomeLikelihood* clone() const
  {
    return new JointTraitChromosomeLikelihood(*this);
  }

  JointTraitChromosomeLikelihood(const JointTraitChromosomeLikelihood& sd);

  //void optimizeTraitModel(double tol, uint numOfIterations, std::vector<string> &traitParamNames, std::vector<string> &fixedParamsTrait);
  //void optimizeChromosomeModel(DerivableSecondOrder* f, double tol, uint numOfIterations, BrentOneDimension* optimizer, std::vector<string> &chromosomeParamNames, std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>> &paramNameAndType);
  void optimizeChromosomeModel(double tol, uint numOfIterations, std::vector<string> &chromosomeParamNames, std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>> &paramNameAndType, std::map<uint, uint> &baseNumberBounds, const string &baseNumOptimizationMethod);
  void optimize(double tol, uint numOfIterations, uint numberOfIterationsPerOneOptimization, std::map<uint, uint> &baseNumberBounds, const string &baseNumOptimizationMethod);
  void setSharedParams(std::map<int, std::vector<std::pair<uint, int>>> &sharedParamsChromosome){
    sharedParamsChromosome_ = sharedParamsChromosome;
  }
  void setFixedParams(std::map<uint, vector<int>> &fixedParamsChromosome){
    fixedParamsChromosome_ = fixedParamsChromosome;
  }
  
  
  
        //   static SingleProcessPhyloLikelihood* setHeterogeneousModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, std::map<int, vector<std::pair<uint,int>>>* updatedSharedParams);
        // static SingleProcessPhyloLikelihood* setRandomHeterogeneousModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, double parsimonyBound, std::map<uint, vector<int>> &fixedParams, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams);






  // Should probably write this function somewhere else
  void printLikParameters(unsigned int optimized, vector<string> paramsNames) const;


};
} // end of namespace bpp.
#endif // CHROMEVOL_JOINTTRAITCHROMOSOMELIKELIHOOD_H
