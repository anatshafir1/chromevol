//
// File: JointTraitChromosomeLikelihood.cpp
// Authors:
//   Anat Shafir
// Created: 2023-04-17 11:57:00
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

#include "BaseNumberOptimizer.h"

using namespace std;
using namespace bpp;


BaseNumberOptimizer::BaseNumberOptimizer(std::shared_ptr<LikelihoodCalculationSingleProcess> lik, bool optimizeBaseNumber, string baseNumOptimizationMethod):
    seqData_(dynamic_cast<const VectorSiteContainer*>(lik->getData())),
    optimizeBaseNumber_(optimizeBaseNumber),
    baseNumOptimizationMethod_(baseNumOptimizationMethod),
    baseNumberUpperBound_() {}

BaseNumberOptimizer::BaseNumberOptimizer(std::shared_ptr<LikelihoodCalculationSingleProcess> lik, bool optimizeBaseNumber, string baseNumOptimizationMethod, std::map<uint, uint> &baseNumberUpperBounds):
    seqData_(dynamic_cast<const VectorSiteContainer*>(lik->getData())),
    optimizeBaseNumber_(optimizeBaseNumber),
    baseNumOptimizationMethod_(baseNumOptimizationMethod),
    baseNumberUpperBound_(baseNumberUpperBounds) {}

BaseNumberOptimizer::BaseNumberOptimizer(const VectorSiteContainer* vsc, std::map<uint, uint> &baseNumberUpperBounds){
    optimizeBaseNumber_ = false;
    seqData_= vsc;
    baseNumberUpperBound_ = baseNumberUpperBounds;

}

void BaseNumberOptimizer::optimizeBaseNum(Function* func, vector<string> paramNames, size_t index, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, double lowerBound, 
                                                double upperBound, const string &paramName, ParameterList& params, uint model, std::map<uint, uint> &baseNumberUpperBounds){

    ParameterList updatedSubstitutionParams = params.createSubList(paramNames);
    size_t best_i = (size_t)(params.getParameter(paramName).getValue());
    double f_value = *currentLikelihood;
    
    for (size_t i = 0; i < baseNumCandidates.size(); i++){
        unsigned int baseNum = baseNumCandidates[i];
        if (baseNum > baseNumberUpperBounds[model]){
            break;
        }
        params.getParameter(paramName).setValue((double)baseNum);
        double f_i = func->f(params);
        if (f_i < f_value){
            best_i = baseNum;
            f_value = f_i;
        }
    }
    params.getParameter(paramName).setValue((double)best_i);
    updatedSubstitutionParams.getParameter(paramName).setValue((double)best_i);
    func->f(updatedSubstitutionParams);
    *currentLikelihood = f_value;

}


/**************************************************************************************************************************************/

void BaseNumberOptimizer::setBaseNumberBounds(std::shared_ptr<LikelihoodCalculationSingleProcess> lik, std::map<uint, uint> &baseNumberBounds,  uint* numberOfModelsPtr){
    
    uint numOfModels;
    if (numberOfModelsPtr){
        numOfModels = *numberOfModelsPtr;
    }else{
        numOfModels = static_cast<uint>(lik->getSubstitutionProcess().getNumberOfModels());
    } 
    for (size_t m = 1; m <= numOfModels; m ++){
        auto branchProcess = lik->getSubstitutionProcess().getModel(m);
        baseNumberBounds[static_cast<uint>(m)] = std::dynamic_pointer_cast<const ChromosomeSubstitutionModel>(branchProcess)->getMaxChrRange();
    }
  }
  

/**************************************************************************************************************************************
 * 
*/
void BaseNumberOptimizer::getBaseNumCandidates(vector <unsigned int> &baseNumCandidates, std::map<uint, uint> &baseNumberUpperBounds) const{
    if ((baseNumOptimizationMethod_ != "Brent") && (optimizeBaseNumber_)){
        uint maxBaseNumCandidate = getMaxBaseNumAmongModels(baseNumberUpperBounds);
        fillVectorOfBaseNumCandidates(baseNumCandidates, lowerLimitBaseNumber, maxBaseNumCandidate);

    }

}

/**************************************************************************************************************************************/
uint BaseNumberOptimizer::getMaxBaseNumAmongModels(std::map<uint, uint> baseNumberUpperBound) const{
    auto it = baseNumberUpperBound.begin();
    uint maxBaseNumBound = 0;
    while(it != baseNumberUpperBound.end()){
        if (baseNumberUpperBound[it->first] > maxBaseNumBound){
            maxBaseNumBound = baseNumberUpperBound[it->first];
        }
        it ++;
    }
    return maxBaseNumBound;
}


// /***********************************************************************************/
void BaseNumberOptimizer::fillVectorOfBaseNumCandidates(std::vector <unsigned int> &baseNumCandidates, unsigned int lowerBound, unsigned int upperBound) const{
    if (baseNumOptimizationMethod_ == "Ranges"){
        getAllPossibleChrRanges(baseNumCandidates);

    }
    else if ((baseNumOptimizationMethod_ == "Sequential") || (baseNumCandidates.size() == 0)){

        for (unsigned int chr = (unsigned int)lowerBound; chr <= upperBound; chr++){
            baseNumCandidates.push_back(chr);
        }

    }

}
// /***************************************************************************************/
void BaseNumberOptimizer::getAllPossibleChrRanges(std::vector <unsigned int> &baseNumCandidates) const{
    size_t numOfSequences = seqData_->getNumberOfSequences();
    unsigned int minRange = 0;
    vector <string> sequenceNames = seqData_->getSequenceNames();
    for (size_t i = 0; i < numOfSequences; i++){
        if (i == numOfSequences-1){
            continue;
        }
        BasicSequence seq1 = seqData_->getSequence(sequenceNames[i]);
        int chrNum1 = seq1.getValue(0);
        if (chrNum1 == -1){
            continue;
        }
        for (size_t j = i + 1; j < numOfSequences; j++){
            BasicSequence seq2 = seqData_->getSequence(sequenceNames[j]);
            int chrNum2 = seq2.getValue(0);
            if (chrNum2 == -1){
                continue;
            }
            unsigned int chrRange = (unsigned int)(std::abs(chrNum1 - chrNum2));
            if (chrRange < lowerLimitBaseNumber){
                continue;
            }

            if (!std::count(baseNumCandidates.begin(), baseNumCandidates.end(), chrRange)){
                if (minRange == 0){
                    minRange = chrRange;
                }else{
                    if (chrRange < minRange){
                        minRange = chrRange;
                    }
                }
                baseNumCandidates.push_back(chrRange);

            }

        }
    }
    if (minRange > lowerLimitBaseNumber){
        for (unsigned int i = lowerLimitBaseNumber; i < minRange; i++){
            baseNumCandidates.push_back(i);
        }

    }

}
