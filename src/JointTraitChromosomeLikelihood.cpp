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

#include "JointTraitChromosomeLikelihood.h"

using namespace std;
using namespace bpp;

JointTraitChromosomeLikelihood::JointTraitChromosomeLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool expectedHistory, bool weightedFrequencies, size_t numOfMappings, std::string traitModel, std::vector<int> &rateChangeType, VectorSiteContainer* chromosomeVsc, std::vector<unsigned int> &baseNumberCandidates, bool ml, bool inCollection):
    AbstractPhyloLikelihood(context),
    JointPhyloLikelihood(context, pC, expectedHistory, weightedFrequencies, numOfMappings, ml, inCollection),
    //BaseNumberOptimizer(std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(tempLik_->getLikelihoodCalculation()), optimizeBaseNumber, baseNumOptimizationMethod, baseNumberUpperBound),
    
    //Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool expectedHistory, bool weightedFrequencies, size_t numOfMappings, bool inCollection = true);
    traitModel_(traitModel),
    fixedParamsChromosome_(),
    fixedParamsTrait_(),
    chromosomeRateChangeType_(rateChangeType),
    baseNumberCandidates_(baseNumberCandidates) ,
    chromosomeVsc_(chromosomeVsc),
    sharedParamsChromosome_()
    {
      // auto lik = std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(tempLik_->getLikelihoodCalculation());
      // if (optimizeBaseNumber){
      //   setBaseNumberBounds(lik, baseNumberUpperBound_);
      //   // getBaseNumCandidates(baseNumberCandidates_, baseNumberUpperBound_);

      // }

    }


JointTraitChromosomeLikelihood::JointTraitChromosomeLikelihood(const JointTraitChromosomeLikelihood& sd):
    AbstractPhyloLikelihood(sd),
    JointPhyloLikelihood(sd),
    traitModel_(sd.traitModel_),
    fixedParamsChromosome_(sd.fixedParamsChromosome_),
    fixedParamsTrait_(sd.fixedParamsTrait_),
    chromosomeRateChangeType_(sd.chromosomeRateChangeType_),
    baseNumberCandidates_(sd.baseNumberCandidates_),
    chromosomeVsc_(sd.chromosomeVsc_),
    sharedParamsChromosome_(sd.sharedParamsChromosome_) {}


// this function can be abstract here, and then redefined by the derived classes
//  void JointTraitChromosomeLikelihood::optimizeTraitModel(double tol, uint numOfIterations, std::vector<string> &traitParamNames, std::vector<string> &fixedParamsTrait){
//   DerivableSecondOrder* f = this;
//   BrentOneDimension* optimizer = new BrentOneDimension(f);
//   optimizer->setVerbose(1);
//   optimizer->setProfiler(0);
//   optimizer->setMessageHandler(0);
//   optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
//   optimizer->setMaximumNumberOfEvaluations(100);
//   // setting bracketing for Brent optimization
//   optimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);
//   traitOptimization_ = true;
//   ParameterList params;
//   double currentLikelihood = getValue();
//   double prevLikelihood;
//   std::cout << "Optimizing trait model ...." << std::endl;
//   for (size_t i = 0; i < numOfIterations; i++){
//     std::cout << "Iteration #" << i << std::endl;
//     for (size_t j = 0; j < traitParamNames.size(); j++){
//         prevLikelihood = currentLikelihood;
//         params = getParameters();
//         const string nameOfParam = traitParamNames[j];
//         std::cout << "Previous value of " << nameOfParam  << " is: "+ std::to_string(params.getParameter(nameOfParam).getValue()) << std::endl;
//         if (std::count(fixedParamsTrait.begin(), fixedParamsTrait.end(), nameOfParam)){
//             continue;
//         }
//         Parameter param = params.getParameter(nameOfParam);
//         std::cout << "Parameter name is: " << nameOfParam << std::endl;
//         std::shared_ptr<IntervalConstraint> bounds = std::dynamic_pointer_cast<IntervalConstraint>(getParameter(nameOfParam).getConstraint());
                 
//         if ((i == 1) & (numOfIterations > 2)){
//             optimizer->getStopCondition()->setTolerance(tol* 2);
//         }else{
//             optimizer->getStopCondition()->setTolerance(tol);
//         }
//         optimizer->setInitialInterval(bounds->getLowerBound()+0.00001, bounds->getUpperBound()-0.00001); 
           
//         optimizer->init(params.createSubList(param.getName()));
//         currentLikelihood = optimizer->optimize();
//         std::cout << "parameter value after optimization "+ std::to_string(getParameter(param.getName()).getValue()) << std::endl;
//         std::cout << "***" << std::endl;

//     }
//     printLikParameters(1, traitParamNames);
//     if (std::abs(prevLikelihood-currentLikelihood) < tol){
//      break;
//     }
//   }
//   delete optimizer;

// }
void JointTraitChromosomeLikelihood::optimizeChromosomeModel(double tol, uint numOfIterations, std::vector<string> &chromosomeParamNames, std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>> &paramNameAndType, std::map<uint, uint> &baseNumberBounds, const string &baseNumberOptimizationMethod){
    //template <class T> unsigned int ChromosomeNumberOptimizer::optimizeModelParameters(T* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, std::vector<string>* chromosomeParamNames)
    ChromosomeNumberOptimizer* chrOptimizer = new ChromosomeNumberOptimizer(static_cast<JointPhyloLikelihood*>(this),
                dynamic_cast<const ChromosomeAlphabet*>(tempLik_->getData()->getAlphabet()),chromosomeVsc_, baseNumberBounds, sharedParamsChromosome_, fixedParamsChromosome_);
    vector<uint> numOfPoints = {1};
    vector<uint> numOfItarationsVec = {numOfIterations};
    chrOptimizer->initOptimizer(numOfPoints, numOfItarationsVec, numOfPoints, numOfItarationsVec, ChromEvolOptions::optimizationMethod_, baseNumberOptimizationMethod,
        tol, ChromEvolOptions::standardOptimization_, ChromEvolOptions::BrentBracketing_, 
        ChromEvolOptions::probsForMixedOptimization_);
    //string* textToPrint = 0;
    //setBaseNumberBounds(std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(tempLik_->getLikelihoodCalculation()), baseNumberUpperBound_);
    //getBaseNumCandidates(baseNumberCandidates_, baseNumberUpperBound_);
    JointPhyloLikelihood* jointLik = static_cast<JointPhyloLikelihood*>(this);
    chrOptimizer->optimizeModelParameters<JointPhyloLikelihood>(jointLik, tol, numOfIterations, baseNumberCandidates_, &(ChromEvolOptions::sharedParameters_), &(ChromEvolOptions::fixedParams_), 0, baseNumberBounds, &chromosomeParamNames);
    delete chrOptimizer;
    

}

// void JointTraitChromosomeLikelihood::optimizeChromosomeModel(DerivableSecondOrder* f, double tol, uint numOfIterations, BrentOneDimension* optimizer, std::vector<string> &chromosomeParamNames, std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>> &paramNameAndType){
//   traitOptimization_ = false;
//   ParameterList params;
//   Function* func = this;
//   double currentLikelihood = getValue();
//   double prevLikelihood;
//   size_t startCompositeParams = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();
//   // starting iterations of optimization
//   std::cout << "Optimizing chromosome number model ...." << std::endl;
//   for (size_t i = 0; i < numOfIterations; i++){
//     std::cout << "Iteration #" << i << endl;
//     for (size_t j = 0; j < chromosomeParamNames.size(); j ++){
//         prevLikelihood = currentLikelihood;
//         params = getParameters();
//         const string nameOfParam = chromosomeParamNames[j];
//         std::cout << "Previous value of " << nameOfParam  << " is: "+ std::to_string(params.getParameter(nameOfParam).getValue()) << std::endl;
//         int rateParamType = paramNameAndType[nameOfParam].first;
            
//         if (std::count(fixedParamsChromosome_[paramNameAndType[nameOfParam].second].begin(), fixedParamsChromosome_[paramNameAndType[nameOfParam].second].end(), rateParamType)){
//             continue;
//         }
//         //int rateCompositeParamType;
//         double lowerBound;
//         double upperBound;
//         // param names corresponding to the parameter type
//         std::vector<string> paramsNames = typeWithParamNames[rateParamType][paramNameAndType[nameOfParam].second];
//         Parameter param = params.getParameter(nameOfParam);
//         auto it = std::find(paramsNames.begin(), paramsNames.end(), nameOfParam);
//         if (it == paramsNames.end()){
//             throw Exception("JointTraitChromosomeLikelihood::optimizeModelParametersOneDimension(): index out of range!");
//         }
//         size_t index = it - paramsNames.begin();
//         if (rateParamType != static_cast<int>(ChromosomeSubstitutionModel::BASENUM)){
//             ChromosomeNumberDependencyFunction::FunctionType funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(chromosomeRateChangeType_[rateParamType-startCompositeParams]);
//             ChromosomeNumberDependencyFunction* functionOp = compositeParameter::setDependencyFunction(funcType);
//             auto alphabet = dynamic_cast<const ChromosomeAlphabet*>(std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(tempLik_->getLikelihoodCalculation())->getData()->getAlphabet());
//             functionOp->setDomainsIfNeeded(alphabet->getMin(), alphabet->getMax());
//             functionOp->updateBounds(params, paramsNames, index, &lowerBound, &upperBound, alphabet->getMax());
//             functionOp->updateBounds(f, nameOfParam, lowerBound, upperBound);
//             delete functionOp;
//             std::shared_ptr<IntervalConstraint> intervalFuncUpdated = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(nameOfParam).getConstraint());
//             std::shared_ptr<IntervalConstraint> intervalFuncUpdatedTL = dynamic_pointer_cast<IntervalConstraint>(getParameter(nameOfParam).getConstraint());

//         }else{
//             // baseNumber parameter
//             if (!std::count(fixedParamsChromosome_[paramNameAndType[nameOfParam].second].begin(), fixedParamsChromosome_[paramNameAndType[nameOfParam].second].end(), ChromosomeSubstitutionModel::BASENUM)){
//                 optimizeBaseNum(func, chromosomeParamNames,j, baseNumberCandidates_, &currentLikelihood, lowerBound, upperBound, nameOfParam, params, paramNameAndType[nameOfParam].second, baseNumberUpperBound_);
//                 std::cout << "parameter value after optimization "+ std::to_string(getParameter(param.getName()).getValue()) << std::endl;
//                 continue;
                
//             }
//             //optimizeBaseNum(Function* func, vector<string> paramNames, size_t index, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, double lowerBound, 
//                                                 //double upperBound, const string &paramName, ParameterList& params, uint model, std::map<uint, uint> &baseNumberUpperBounds)
//         }
//         std::cout << "Parameter name is: " << nameOfParam << std::endl;
                 
//         if ((i == 1) & (numOfIterations > 2)){
//             optimizer->getStopCondition()->setTolerance(tol* 2);
//         }else{
//             optimizer->getStopCondition()->setTolerance(tol);
//         }
//         if (rateParamType != static_cast<int>(ChromosomeSubstitutionModel::BASENUM)){
//             optimizer->setInitialInterval(lowerBound + 1e-10, upperBound);
//         }else{
//             optimizer->setInitialInterval(lowerBound, upperBound);
//         }            
//         optimizer->init(params.createSubList(param.getName()));
//         currentLikelihood = optimizer->optimize();
//         std::cout << "parameter value after optimization "+ std::to_string(getParameter(param.getName()).getValue()) << std::endl;
//         std::cout << "***" << std::endl;
//     }
//     printLikParameters(1, chromosomeParamNames);
        
//     if (std::abs(prevLikelihood-currentLikelihood) < tol){
//         break;
//     }
//   }
       
// }

void JointTraitChromosomeLikelihood::printLikParameters(unsigned int optimized, vector<string> paramsNames) const{

    if (optimized == 0){
        std::cout << "Initial likelihood is : "<< getValue() << std::endl;
        std::cout << "Initial trait likelihood is: " << getAbstractPhyloLikelihood(getNumbersOfPhyloLikelihoods()[0])->getValue() << std::endl;
        std::cout << "Initial chromosome likelihood is: " << tempLik_->getValue() << std::endl;
        std::cout << "Parameters are:" << std::endl;
    }else{
        std::cout << "Optimized likelihood is : " << getValue() << std::endl;
        std::cout << "Optimized trait likelihood is: " << getAbstractPhyloLikelihood(getNumbersOfPhyloLikelihoods()[0])->getValue() << std::endl;
        std::cout << "Optimized chromosome likelihood is: " << tempLik_->getValue() << std::endl;
        std::cout << "Optimized parameters are:"<< std::endl;
    }   
    for (int i = 0; i < (int)(paramsNames.size()); i++){
        if (paramsNames[i].find("Chromosome.baseNum_") != std::string::npos){
           std::cout << paramsNames[i] << " = " << (int)(getParameter(paramsNames[i]).getValue()) << std::endl;
        }else{
            std::cout << paramsNames[i] << " = " << getParameter(paramsNames[i]).getValue() << std::endl;

        }    
    }
    std::cout <<  "***" << std::endl;

}

  



void JointTraitChromosomeLikelihood::optimize(double tol, uint numOfIterations, uint numberOfIterationsPerOneOptimization, std::map<uint, uint>& baseNumberBounds, const string &baseNumberOptimizationMethod){

    // initializing the likelihood values
    
    // separating between trait parameters and chromosome number evolution parameters
    std::map<uint, std::vector<std::string>> paramsPerModel;
    LikelihoodUtils::separateBetweenModels(this, traitModel_, paramsPerModel);
    std::vector<string> traitParamNames = paramsPerModel[1];
    std::vector<string> chromosomeParamNames = paramsPerModel[2];
    // setting maps of parameter type and the corresponding parameters, and vice versa
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model

    LikelihoodUtils::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, chromosomeParamNames, 0);
    for (size_t i = 0; i< numOfIterations; i++){
        std::cout << "Iteration #" << i << std::endl;
        traitOptimization_ =false;
        //std::cout << "Cycle " << i << std::endl;
        std::cout << "Joint likelihood value (before trait) " << getValue() << std::endl;
        std::cout << "trait likelihood value (before trait) " << getAbstractPhyloLikelihood(nPhylo_[0])->getValue() << std::endl;
        std::cout << "Chromosome likelihood value (before trait) " << tempLik_->getValue() << std::endl;
        std::cout << "*** Optimizing trait parameters ... " << std::endl;

        traitOptimization_ =true;
        optimizeTraitModel<JointPhyloLikelihood>(static_cast<JointPhyloLikelihood*>(this), tol, numberOfIterationsPerOneOptimization, traitParamNames, fixedParamsTrait_);
        

        //optimizeChromosomeModel(f, tol, numberOfIterationsPerOneOptimization, optimizer, chromosomeParamNames, typeWithParamNames, paramNameAndType);
        traitOptimization_ = false;
        std::cout << "Joint likelihood value (after trait) " << getValue() << std::endl;
        std::cout << "trait likelihood value (after trait)" << getAbstractPhyloLikelihood(nPhylo_[0])->getValue() << std::endl;
        std::cout << "Chromosome likelihood value (after trait)" << tempLik_->getValue() << std::endl;
        std::cout << "*** Optimizing chromosome parameters ... " << std::endl;

        optimizeChromosomeModel(tol, numberOfIterationsPerOneOptimization, chromosomeParamNames, typeWithParamNames, paramNameAndType, baseNumberBounds, baseNumberOptimizationMethod);
        std::cout << "Joint likelihood value (after chromosome) " << getValue() << std::endl;
        std::cout << "trait likelihood value (after chromosome)" << getAbstractPhyloLikelihood(nPhylo_[0])->getValue() << std::endl;
        std::cout << "Chromosome likelihood value (after chromosome)" << tempLik_->getValue() << std::endl;
        std::cout << "** End of Iteration **" << std::endl;

    }
    


}

    // std::map<uint, uint> maxBaseNumTransition = (ChromEvolOptions::simulateData_ || ChromEvolOptions::useMaxBaseTransitonNumForOpt_) ? ChromEvolOptions::maxBaseNumTransition_ : chrRange_;
    // ChromosomeNumberOptimizer* opt = new ChromosomeNumberOptimizer(tree_, alphabet_, vsc_, maxBaseNumTransition);
    // //initialize all the optimization specific parameters
    // opt->initOptimizer(ChromEvolOptions::OptPointsNum_, ChromEvolOptions::OptIterNum_, ChromEvolOptions::OptPointsNumNextRounds_, ChromEvolOptions::OptIterNumNextRounds_, ChromEvolOptions::optimizationMethod_, ChromEvolOptions::baseNumOptimizationMethod_,
    //     ChromEvolOptions::tolerance_, ChromEvolOptions::standardOptimization_, ChromEvolOptions::BrentBracketing_, 
    //     ChromEvolOptions::probsForMixedOptimization_);
    // //optimize models
    // time_t t1;
    // time(&t1);
    // time_t t2;
    // std::map<uint, vector<uint>> mapOfNodesWithRoot = ChromEvolOptions::mapModelNodesIds_;
    // mapOfNodesWithRoot[1].push_back(tree_->getRootIndex());
    // auto modelAndRepresentitives = LikelihoodUtils::findMRCAForEachModelNodes(tree_, mapOfNodesWithRoot);
    // opt->setInitialModelRepresentitives(modelAndRepresentitives);
    // if (ChromEvolOptions::parallelization_){
    //     opt->optimizeInParallel(complexParamsValues, parsimonyBound, ChromEvolOptions::rateChangeType_, ChromEvolOptions::seed_, ChromEvolOptions::OptPointsNum_[0], ChromEvolOptions::fixedFrequenciesFilePath_, ChromEvolOptions::fixedParams_ ,ChromEvolOptions::mapModelNodesIds_);

    // }else{
    //     opt->optimize(complexParamsValues, parsimonyBound, ChromEvolOptions::rateChangeType_, ChromEvolOptions::seed_, ChromEvolOptions::OptPointsNum_[0], ChromEvolOptions::fixedFrequenciesFilePath_, ChromEvolOptions::fixedParams_ ,ChromEvolOptions::mapModelNodesIds_);

    // }