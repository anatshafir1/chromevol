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

#ifndef CHROMEVOL_BASENUMBEROPTIMIZER_H
#define CHROMEVOL_BASENUMBEROPTIMIZER_H


#include <Bpp/Numeric/AutoParameter.h>

//from bpp-seq
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>


//from bpp-phyl

#include "ChromosomeSubstitutionModel.h"
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
// From Seqlib:
#include <vector>
#include <map>
#include <utility>
#include <string>
namespace bpp
{
/**
 * @brief The BaseNumberOptimizer class: is needed for optimizing base number in chromosome number models
 *
 *
 */ 

class BaseNumberOptimizer
{
  protected:
    const VectorSiteContainer* seqData_;
    bool optimizeBaseNumber_;
    string baseNumOptimizationMethod_;
    mutable std::map<uint, uint> baseNumberUpperBound_;



public:
  BaseNumberOptimizer(std::shared_ptr<LikelihoodCalculationSingleProcess> lik, bool optimizeBaseNumber, string baseNumOptimizationMethod);
  BaseNumberOptimizer(std::shared_ptr<LikelihoodCalculationSingleProcess> lik, bool optimizeBaseNumber, string baseNumOptimizationMethod, std::map<uint, uint> &baseNumberUpperBounds);
  BaseNumberOptimizer(const VectorSiteContainer* vsc, std::map<uint, uint> &baseNumberUpperBounds);

  BaseNumberOptimizer(const BaseNumberOptimizer &bn):
  seqData_(bn.seqData_),
  optimizeBaseNumber_(bn.optimizeBaseNumber_),
  baseNumOptimizationMethod_(bn.baseNumOptimizationMethod_),
  baseNumberUpperBound_(bn.baseNumberUpperBound_)
  {}

  ~BaseNumberOptimizer() {}

  BaseNumberOptimizer* clone() const
  {
    return new BaseNumberOptimizer(*this);
  }
  void setBaseNumberOptimizationMethod(string baseNumOptimizationMethod){
    baseNumOptimizationMethod_ = baseNumOptimizationMethod;
  }
  void setIfBaseNumberOptimized(bool optimizeBaseNumber){
    optimizeBaseNumber_ = optimizeBaseNumber;
  }



protected:
  void optimizeBaseNum(Function* func, vector<string> paramNames, size_t index, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, double lowerBound, 
                                                double upperBound, const string &paramName, ParameterList& params, uint model, std::map<uint, uint> &baseNumberUpperBounds);
  // will be invoked separately (not in the constructor, because it might be also given as argument)
  void setBaseNumberBounds(std::shared_ptr<LikelihoodCalculationSingleProcess> lik, std::map<uint, uint> &baseNumberBounds, uint* numOfModelsPtr = 0);
  void getBaseNumCandidates(vector <unsigned int> &baseNumCandidates, std::map<uint, uint> &baseNumberUpperBounds) const;
  uint getMaxBaseNumAmongModels(std::map<uint, uint> baseNumberUpperBound) const;
  void fillVectorOfBaseNumCandidates(std::vector <unsigned int> &baseNumCandidates, unsigned int lowerBound, unsigned int upperBound) const;
  void getAllPossibleChrRanges(std::vector <unsigned int> &baseNumCandidates) const;



};
} // end of namespace bpp.
#endif // CHROMEVOL_BASENUMBEROPTIMIZER_H
