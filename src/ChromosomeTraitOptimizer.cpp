#include "ChromosomeTraitOptimizer.h"
using namespace bpp;

std::map <uint, std::vector<uint>> ChromosomeTraitOptimizer::getNodesForEachModel(std::shared_ptr<PhyloTree> expectedMapping, StochasticMapping* stm){
  std::map <uint, std::vector<uint>> nodeModels;
  auto nodes = expectedMapping->getAllNodes();
  for (size_t i = 0; i < nodes.size(); i ++){
    uint nodeId = expectedMapping->getNodeIndex(nodes[i]);
    int nodeState = stm->getNodeState(nodes[i]);
    if (!(expectedMapping->isLeaf(nodeId))){
      auto sons = expectedMapping->getSons(nodeId);
      for (size_t j = 0; j < sons.size(); j++){
        // add +1 so that it will work for the general function that constructs the heterogeneous model
        nodeModels[static_cast<uint>(nodeState)+1].push_back(sons[j]);
      }
    }
  }
  return nodeModels;
}


vector<double> ChromosomeTraitOptimizer::getFrequenciesFromMapOfParams(std::map<string, double> &traitModelParams, bool random){
  vector<double> freqs;
  double sumOfFreqs = 0;
  for (size_t i = 0; i < static_cast<size_t>(ChromEvolOptions::numberOfTraitStates_); i++){
    if (random){
      if (i == static_cast<size_t>(ChromEvolOptions::numberOfTraitStates_)-1){
        freqs.push_back(1-sumOfFreqs);
      }else{
        auto freq = RandomTools::giveRandomNumberBetweenTwoPoints(0.05, std::min(0.95, 1-sumOfFreqs));
        sumOfFreqs += freq;
        freqs.push_back(freq);
      }
    }else{
      freqs.push_back(traitModelParams["pi"+ std::to_string(i)]);
    }
  }
  return freqs;
}




/*
The initiation of the likelihood function is either with some fixed values or with random values. If the values are random,
it makes sense to start from assigning the same model for both states.
*/

void ChromosomeTraitOptimizer::initJointLikelihood(std::map<string, double> traitModelParams, std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, const PhyloTree* tree, std::map<uint, uint> baseNumberUpperBound, std::vector<double>* rootFreqsTrait, bool random, double maxParsimonyFactor, bool ml){
  std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT;
  std::shared_ptr<NonHomogeneousSubstitutionProcess> nsubProT;
  std::shared_ptr<LikelihoodCalculationSingleProcess> likT;
  //bool traitWeightedFrequencies = true;
  std::shared_ptr<BranchModel> traitModel;
  std::vector<uint> modelNodesTrait;
  std::shared_ptr<DiscreteDistribution> rdistTrait = make_shared<GammaDiscreteRateDistribution>(1, 1.0);
  Context contextT = Context();
  const IntegerAlphabet* traitAlpha = dynamic_cast<const IntegerAlphabet*>(vscTrait_->getAlphabet());
  vector<double> freqVals;
  if (fixedTraitRootFreq_){
    freqVals = getFrequenciesFromMapOfParams(traitModelParams, false);
  }else{
    freqVals = getFrequenciesFromMapOfParams(traitModelParams, random); 

  }
  shared_ptr<IntegerFrequencySet> freqs = make_shared<FullIntegerFrequencySet>(traitAlpha, freqVals);
  std::shared_ptr<CharacterSubstitutionModel> model = LikelihoodUtils::setTraitModel(traitAlpha, freqs);
  
  //ratePerPairModel.setParameterValue("rate_1_0", rateVals(1, 0));
  setParametersToNewTraitModel(traitModelParams, model, freqs, random);

  traitModel = static_pointer_cast<BranchModel>(model);
  std::shared_ptr<FrequencySet> rootFrequenciesFixedRoot;
  
  if (fixedRootTraitState_ >= 0){
    vector<double> rootTraitFreqs;
    rootTraitFreqs.resize(ChromEvolOptions::numberOfTraitStates_);
    for (size_t k = 0; k < ChromEvolOptions::numberOfTraitStates_; k++){
      if (k == fixedRootTraitState_){
        rootTraitFreqs[k] = 1;
      }else{
        rootTraitFreqs[k] = 0;
      }
    }
    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(traitModel->getStateMap(), false)), rootTraitFreqs);
    rootFrequenciesFixedRoot = std::shared_ptr<FrequencySet>(rootFreqsFixed->clone());
    subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree, rootFrequenciesFixedRoot);

  }else{
    subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);
  }


  //subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);
  getTraitNodes(modelNodesTrait);
    
  subProT->addModel(std::shared_ptr<CharacterSubstitutionModel>(model->clone()), modelNodesTrait);
  nsubProT = std::shared_ptr<NonHomogeneousSubstitutionProcess>(subProT->clone());
  string prefix  = traitModel->getNamespace();
  //LikelihoodUtils::aliasTraitParams(std::shared_ptr<NonHomogeneousSubstitutionProcess> process, int &numOfTraitConstraints, string &prefix, std::unordered_map<string, string> &sharedTraitParams);
  LikelihoodUtils::aliasTraitParams(nsubProT, ChromEvolOptions::numOfTraitConstraints_, prefix, ChromEvolOptions::sharedTraitParams_);



  Context* context = new Context();
  likT = std::make_shared<LikelihoodCalculationSingleProcess>(contextT, *vscTrait_->clone(), *nsubProT, weightedTraitRootFreqs_);
  auto phyloT = std::make_shared<SingleProcessPhyloLikelihood>(contextT, likT);
  auto logLikVal = phyloT->getValue(); // calculate now because it will be needed for stochastic mapping
  std::cout << "Initial trait likelihood in joint model: " << std::endl;

  // get expected tree for the chromosome model
  StochasticMapping* stm = new StochasticMapping(likT, numberOfStochasticMappings_);
  std::map<uint, std::vector<size_t>> mlAncestors;
  if (ml){
    mlAncestors = JointPhyloLikelihood::getMLAncestralReconstruction(phyloT.get(), phyloT->getParameters(), 2);
    stm->setMLAncestors(&mlAncestors);
  }
  auto treeChr = stm->createExpectedMappingHistory(numberOfStochasticMappings_);
  // creating the chromosome number model
  std::map<uint, vector<uint>> mapModelNodesIds = getNodesForEachModel(treeChr, stm);
  delete stm;
  std::shared_ptr<ParametrizablePhyloTree> parTreeChr;
  std::vector<std::shared_ptr<ChromosomeSubstitutionModel>> models;
  SubstitutionProcess* subProcess;
  if (random){
    subProcess = LikelihoodUtils::setRandomChromosomeSubstitutionModel(treeChr.get(), vscChr_, alphabetChr_, baseNumberUpperBound_, mapModelNodesIds, modelParams, ChromEvolOptions::numberOfTraitStates_, maxParsimonyFactor, fixedParams_, &sharedParams_, weightedFreqs_, &models, parTreeChr);

  }else{
    subProcess = LikelihoodUtils::setChromosomeSubstitutionModel(treeChr.get(), vscChr_, alphabetChr_, baseNumberUpperBound_, mapModelNodesIds, modelParams, ChromEvolOptions::numberOfTraitStates_, &sharedParams_, weightedFreqs_, &models, parTreeChr);

  }
  auto parTree2 =  std::make_shared<ParametrizablePhyloTree>(*treeChr);
  

  auto modelColl=std::make_shared<SubstitutionProcessCollection>();
  
  modelColl->addModel(traitModel, 1);
  for (size_t i = 0; i < models.size(); i++){
    modelColl->addModel(models[i], i+2);

  }

  // if (!traitWeightedFrequencies){
  //   modelColl->addFrequencies(rootFrequenciesTrait, 1);

  // }

  
  //auto root_freq_cloned = rootFreqs2->clone();
  //auto rootFreqs_chr = std::shared_ptr<FrequencySet>(root_freq_cloned);
  //modelColl->addFrequencies(rootFreqs_chr, 2);
  std::shared_ptr<DiscreteDistribution> rdistChr = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
  
  modelColl->addDistribution(rdistTrait, 1);
  modelColl->addDistribution(rdistChr, 2);

  modelColl->addTree(parTree, 1);
  modelColl->addTree(parTree2, 2);
  vector<uint> vP1m1 = modelNodesTrait;

  map<size_t, Vuint> mModBr1;
  mModBr1[1]=vP1m1;
 

  modelColl->addSubstitutionProcess(1, mModBr1, 1, 1);
                                   
  map<size_t, Vuint> mModBr2;

  for (size_t i = 0; i < models.size(); i++){
    mModBr2[i+2]= mapModelNodesIds[i+1];
  }


  modelColl->addSubstitutionProcess(2, mModBr2, 2, 2);

  // Now setting the joint likelihood object
    // Likelihoods
  SubstitutionProcess* sP1c=nsubProT->clone();
  SubstitutionProcess* sP2c= subProcess;//subProcess->clone(); // no need to clone again- it was already cloned.
  //delete subProcess;
  auto pc(std::make_shared<PhyloLikelihoodContainer>(*context, *modelColl));



  auto lik1_j = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *(vscTrait_->clone()), *sP1c, weightedTraitRootFreqs_);

  pc->addPhyloLikelihood(1, new SingleProcessPhyloLikelihood(*context, lik1_j));
    
  auto lik2_j = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *(vscChr_->clone()), *sP2c, weightedFreqs_);

  pc->addPhyloLikelihood(2, new SingleProcessPhyloLikelihood(*context, lik2_j));
  //   JointTraitChromosomeLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool expectedHistory, bool weightedFrequencies, size_t numOfMappings, std::string traitModel, std::vector<int> &rateChangeType, VectorSiteContainer* chromosomeVsc, std::vector<unsigned int> &baseNumberCandidates, bool inCollection = true);
  JointTraitChromosomeLikelihood* jl = new JointTraitChromosomeLikelihood(*context, pc, true, true, numberOfStochasticMappings_, traitModel_, ChromEvolOptions::rateChangeType_, vscChr_->clone(), baseNumberCandidates_, ml);
  jl->setStochasticMappingTree(treeChr);
  jl->setSharedParams(sharedParams_);
  jl->setFixedParams(fixedParams_);
  vectorOfJointLikelohoods_.push_back(jl);
  std::cout << "Initial likelihood is " << jl->getValue() << std::endl;
  auto params = jl->getSubstitutionModelParameters();
  std::cout <<"Parameters are:" << std::endl;
  for (size_t i = 0; i < params.size(); i++){
    std::cout <<"\t" << params[i].getName() << ": " << params[i].getValue() << std::endl;

  }

}
// std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, const PhyloTree* tree, const VectorSiteContainer* vsc_chr, const VectorSiteContainer* vsc_trait, const ChromosomeAlphabet* alphabet_chr, std::vector<double>* traitModelParams, std::map<uint, uint> baseNumberUpperBound, std::map<uint, pair<int, std::map<int, std::vector<double>>>>* modelParams_chr, std::map<int, vector<std::pair<uint,int>>>* updatedSharedParams_chr, std::vector<double>* rootFreqsChr, std::vector<double>* rootFreqsTrait, bool random

void ChromosomeTraitOptimizer::initMultipleLikelihoodPoints(std::map<string, double> &traitModelParams, std::map<uint, std::pair<int, std::map<int, vector<double>>>> &modelParams, const PhyloTree* tree, std::map<uint, uint> baseNumberUpperBound, std::vector<double>* rootFreqsTrait, bool ml){
    std::cout << "******************************************************" << std::endl;
    std::cout << "*  Optimizing joint likelihood function  *" << std::endl;
    std::cout << "******************************************************" << std::endl;
    size_t index = min((int)numOfPoints_.size()-1, 1);
    vectorOfJointLikelohoods_.reserve(numOfPoints_[index]);
    // If base number is one of the parameters
    cout <<"##################################" << endl;
    cout << "*********  cycle 0  **************"<<endl;  
    for (size_t n = 0; n < numOfPoints_[0]; n++){
        std::cout << "Starting cycle with Point #" << n <<"...."<<endl;
        if (n == 0){
          initJointLikelihood(traitModelParams, modelParams, tree, baseNumberUpperBound, rootFreqsTrait, false, 0, ml);
        }else{
          double factor;
          if (n > 1){
            factor = parsimonyBound_ * (1+(0.1*(double)n));

          }else{
            factor = parsimonyBound_ * static_cast<double>(n);
          }
          
          initJointLikelihood(traitModelParams, modelParams, tree, baseNumberUpperBound, rootFreqsTrait, true, factor, ml);
    
        }   
    }

    sort(vectorOfJointLikelohoods_.begin(), vectorOfJointLikelohoods_.end(), compareJointLikValues);



}
void ChromosomeTraitOptimizer::fillTraitParameters(const string &traitModel, vector<string> &parameterNames, vector<string> &traitParamNames) const{
  for (size_t i = 0; i < parameterNames.size(); i++){
    if (parameterNames[i].find(traitModel) != std::string::npos){
      traitParamNames.push_back(parameterNames[i]);
    } 
    
  }

}
/*********************************************************************************************
 * 
 * 
*/
size_t ChromosomeTraitOptimizer::getNumberOfParametersJointLikelihood(){
  std::map<uint, std::vector<std::string>> paramsPerModel;
  LikelihoodUtils::separateBetweenModels(vectorOfJointLikelohoods_[0], traitModel_, paramsPerModel);
  std::vector<string> traitParamNames = paramsPerModel[1];
  std::vector<string> chromosomeParamNames = paramsPerModel[2];
  size_t numOfAllParams = traitParamNames.size() + chromosomeParamNames.size();
  size_t fixedChromosomeParams = LikelihoodUtils::getNumberOfFixedParams(vectorOfJointLikelohoods_[0]->getPhylo2(), fixedParams_);
  size_t fixedRootFreqParams = 0;
  if (fixedTraitRootFreq_){
    fixedRootFreqParams += ChromEvolOptions::numberOfTraitStates_-1;

  }
  size_t numOfActualParams = numOfAllParams - fixedTraitParams_.size() - fixedRootFreqParams - fixedChromosomeParams;
  return numOfActualParams;
}
/*********************************************************************************************/
size_t ChromosomeTraitOptimizer::getNumberOfParametersNull(){
  size_t fixedChromosomeParams = LikelihoodUtils::getNumberOfFixedParams(optimizedChromosomeLikelihood_, fixedParams_);
  auto params = vectorOfLikelihoodsTrait_[0]->getSubstitutionModelParameters();
  vector<string> parameterNames;
  LikelihoodUtils::fixSuffixForJointLikParamNames(params, parameterNames);
  std::vector<string> traitParamNames;
  string traitPrefix = vectorOfLikelihoodsTrait_[0]->getSubstitutionProcess().getModel(1)->getNamespace();
  fillTraitParameters(traitPrefix, parameterNames, traitParamNames);
  size_t totalNumOfTraitParams = traitParamNames.size();
  size_t totalNumberOfChromosomeParameters = optimizedChromosomeLikelihood_->getSubstitutionModelParameters().size();
  size_t fixedRootFreqParams = 0;
  if (fixedTraitRootFreq_){
    fixedRootFreqParams += ChromEvolOptions::numberOfTraitStates_-1;

  }
  size_t numOfActualParams = totalNumOfTraitParams + totalNumberOfChromosomeParameters - fixedChromosomeParams - fixedTraitParams_.size() - fixedRootFreqParams;
  return numOfActualParams;
}
/*********************************************************************************************/

void ChromosomeTraitOptimizer::optimizeIndependentLikelihood(){
  // since we assume that the independent chromosome likelihood model has been already optimized,
  // I can continue with optimizing only trait model.
  // This function should be used after the function initTraitLikelihoods()
  std::cout << "******************************************************" << std::endl;
  std::cout << "*  Optimizing trait independent likelihood function  *" << std::endl;
  std::cout << "******************************************************" << std::endl;
  auto params = vectorOfLikelihoodsTrait_[0]->getSubstitutionModelParameters();
  std::vector<string> parameterNames = params.getParameterNames();
  auto fixed_params_names = getTraitFixedParamFullNames(parameterNames, fixedTraitParams_, fixedTraitRootFreq_);

  for (size_t i = 0; i < numOfIterations_.size(); i++){
    if (i > 0){
      cout <<"##################################" << endl;
      cout << "*********  cycle " << i << "  **************"<<endl;  

    }
    for (size_t j = 0; j < numOfPoints_[i]; j++){
      //If the number of optimization iterations is larger than zero, optimize the number of times as specified
      std::cout << "Optimizing point #" << j << "..." << std::endl;
      if (numOfIterations_[i] > 0){
        // (T* lik, double tol, uint numOfIterations, std::vector<string> &traitParamNames, std::vector<string> &fixedParamsTrait)
        optimizeTraitModel<SingleProcessPhyloLikelihood>(vectorOfLikelihoodsTrait_[j], tolerance_, numOfIterations_[i], parameterNames, fixed_params_names);
        //vectorOfJointLikelohoods_[j]->optimize(tolerance_, 1, numOfIterationInBetween_, baseNumberUpperBound_, baseNumOptimizationMethod_);
      }

    }
    sort(vectorOfLikelihoodsTrait_.begin(), vectorOfLikelihoodsTrait_.end(), LikelihoodUtils::compareLikValues);
    if (i < numOfIterations_.size()-1){
      clearVectorOfLikelihoods(vectorOfLikelihoodsTrait_, numOfPoints_[i+1]);
    }
    for (size_t k = 0; k < vectorOfLikelihoodsTrait_.size(); k++){
      std::cout << "Likelihoods at the end of the cycle are: " << vectorOfLikelihoodsTrait_[k]->getValue() << std::endl;
    }
  }



}
/*******************************************************************************************************/
void ChromosomeTraitOptimizer::initTraitLikelihoods(std::map<string, double> &traitParams){
    size_t index = min((int)numOfPoints_.size()-1, 1);
    vectorOfLikelihoodsTrait_.reserve(numOfPoints_[index]);
    // If base number is one of the parameters
    cout <<"##################################" << endl;
    cout << "*********  cycle 0  **************"<<endl;  
    for (size_t n = 0; n < numOfPoints_[0]; n++){
        std::cout << "Starting cycle with Point #" << n <<"...."<<endl;
        if (n == 0){
          initTraitLikelihood(traitParams, false);
        }else{
          initTraitLikelihood(traitParams, true);
        }   
    }

    sort(vectorOfLikelihoodsTrait_.begin(), vectorOfLikelihoodsTrait_.end(), LikelihoodUtils::compareLikValues);

}
/*********************************************************************************************************/
void ChromosomeTraitOptimizer::initTraitLikelihood(std::map<string, double> &traitParams, bool random){
  std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree_);
  const IntegerAlphabet* traitAlpha = dynamic_cast<const IntegerAlphabet*>(vscTrait_->getAlphabet());
  std::shared_ptr<DiscreteDistribution> rdistTrait = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));

  vector<double> freqVals;
  if (fixedTraitRootFreq_){
    freqVals = getFrequenciesFromMapOfParams(traitParams, false); 

  }else{
    freqVals = getFrequenciesFromMapOfParams(traitParams, random); 

  }
  
  shared_ptr<IntegerFrequencySet> freqs = make_shared<FullIntegerFrequencySet>(traitAlpha, freqVals);
  std::shared_ptr<CharacterSubstitutionModel> characterModel = LikelihoodUtils::setTraitModel(traitAlpha, freqs);
  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT;
  setParametersToNewTraitModel(traitParams, characterModel, freqs, random);
  if (fixedRootTraitState_ >= 0){
    vector<double> rootTraitFreqs;
    rootTraitFreqs.resize(ChromEvolOptions::numberOfTraitStates_);
    for (size_t k = 0; k < ChromEvolOptions::numberOfTraitStates_; k++){
      if (k == fixedRootTraitState_){
        rootTraitFreqs[k] = 1;
      }else{
        rootTraitFreqs[k] = 0;
      }
    }
    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(characterModel->getStateMap(), false)), rootTraitFreqs);
    std::shared_ptr<FrequencySet> rootFrequenciesFixedRoot = std::shared_ptr<FrequencySet>(rootFreqsFixed->clone());
    subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree, rootFrequenciesFixedRoot);

  }else{
    subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);
  }

  
  std::vector<uint> modelNodesTrait;
  getTraitNodes(modelNodesTrait);
  string prefix  = characterModel->getNamespace();
  subProT->addModel(std::shared_ptr<CharacterSubstitutionModel>(characterModel->clone()), modelNodesTrait);
  LikelihoodUtils::aliasTraitParams(subProT, ChromEvolOptions::numOfTraitConstraints_, prefix, ChromEvolOptions::sharedTraitParams_);
  Context* context = new Context();
  auto likT = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vscTrait_->clone(), *(subProT->clone()), weightedTraitRootFreqs_);
  auto phyloT = new SingleProcessPhyloLikelihood(*context, likT);
  vectorOfLikelihoodsTrait_.push_back(phyloT);
  std::cout << "Initial likelihood is: " << phyloT->getValue() << std::endl;
  std::cout << "Model parameters are: " << std::endl;
  auto substitutionModelParams = phyloT->getSubstitutionModelParameters();
  for (size_t i = 0; i < substitutionModelParams.size(); i++){
    std::cout << "\t" << substitutionModelParams[i].getName() << " value is " << substitutionModelParams[i].getValue() << std::endl;
  }


}
/******************************************************************************************/
double ChromosomeTraitOptimizer::calculateFreqs(vector<double> &thetas, size_t &idx) const{
    double res;
    if (idx == 0){
        res = thetas[0];
    }else{
        if (idx < ChromEvolOptions::numberOfTraitStates_-1){
            res = thetas[idx];
            for (size_t i = 0; i < idx; i++){
                res *= (1-thetas[i]);
            }
        }
        else{
            res = 1;
            for (size_t i = 0; i < idx; i++){
                res *= (1-thetas[i]);
            }
        }
    }
    return res;

}
/******************************************************************************************/
std::map<string, double> ChromosomeTraitOptimizer::getTraitMLParamsIndependentLik(){
  auto paramsNames = vectorOfLikelihoodsTrait_[0]->getSubstitutionModelParameters().getParameterNames();
  auto allParams = vectorOfLikelihoodsTrait_[0]->getParameters();
  auto characterModel = (vectorOfLikelihoodsTrait_[0]->getSubstitutionProcess()).getModel(1);
  string modelPrefix = characterModel->getNamespace();
  string suffix = "_1"; // this is always the suffix, because in the trait model we only have one model
  std::map<string, double> traitParams;

  // setting root frequencies
  ValueRef <Eigen::RowVectorXd> rootFreqs = vectorOfLikelihoodsTrait_[0]->getLikelihoodCalculationSingleProcess()->getRootFreqs();
  auto rootFreqsValues =  rootFreqs->getTargetValue();
  Vdouble rootFreqsBpp;
  copyEigenToBpp(rootFreqsValues, rootFreqsBpp);
  size_t numberOfStates = rootFreqsBpp.size();
  if ((!weightedTraitRootFreqs_) && (fixedRootTraitState_ < 0)){
    for (size_t i = 0; i < numberOfStates; i++){
      traitParams["pi"+ std::to_string(i)] = rootFreqsBpp[i];
    }

  }

  vector<double> thetas;
  thetas.resize(ChromEvolOptions::numberOfTraitStates_-1);
  
  for (size_t i = 0; i < paramsNames.size(); i++){
    auto nameWithoutPrefix = (paramsNames[i]).substr(modelPrefix.length());
    auto name = nameWithoutPrefix.substr(0, nameWithoutPrefix.size() - suffix.length());
    auto theta_index = getThetaIndexIfTheta(name);
    if (theta_index >= 0){
      thetas[(size_t)theta_index] = allParams.getParameter(paramsNames[i]).getValue();
    }else{
      traitParams[name] = allParams.getParameter(paramsNames[i]).getValue();
    }   
  }
  if ((weightedTraitRootFreqs_) || (fixedRootTraitState_ >= 0)){
    for (size_t i = 0; i < ChromEvolOptions::numberOfTraitStates_; i++){
      traitParams["pi"+ std::to_string(i)] = calculateFreqs(thetas, i);
    }
  }
  return traitParams;
}

/******************************************************************************************/
int ChromosomeTraitOptimizer::getThetaIndexIfTheta(const std::string& paramName) {
  // Define the regex pattern to match 'theta' followed by a digit
  std::regex pattern(R"(theta(\d))");
  std::smatch match;
  int res = -1;
  // Search for the pattern in the input string
  if (std::regex_search(paramName, match, pattern)) {
    // If a match is found, convert the captured group (the digit) to an integer and return it
    res = std::stoi(match[1].str()) - 1;
  } 
  return res;
}

/******************************************************************************************/
std::map<uint, std::pair<int, std::map<int, vector<double>>>> ChromosomeTraitOptimizer::getChromosomeMLParamsIndependent(uint numOfRequiredModels){
  std::map<uint, std::pair<int, std::map<int, vector<double>>>> chromEvolParams;
  auto &lik = optimizedChromosomeLikelihood_;
  std::map<std::string, int> typeGeneralName;
  LikelihoodUtils::updateWithTypeAndCorrespondingName(typeGeneralName);
  auto allParams = lik->getParameters();
  auto paramNames = lik->getSubstitutionModelParameters().getParameterNames();
  uint numOfModels = static_cast<uint>(lik->getSubstitutionProcess().getNumberOfModels());
  std::map<uint, vector<string>> modelsToNames;
  for (auto &paramName : paramNames){
    for (uint i = 1; i <= numOfModels; i++){
      std::string suffix = "_" + std::to_string(i);
      if (paramName.compare(paramName.length() - suffix.length(), suffix.length(), suffix) == 0){
        modelsToNames[i].push_back(paramName);
        break;
      }

    }
  }

  //////////////////////////////////////
  for (uint i = 1; i <= numOfModels; i++){
    auto &paramsPerModel = modelsToNames[i];
    chromEvolParams[i];
    for (auto &fullName : paramsPerModel){
      auto it = typeGeneralName.begin();
      while (it != typeGeneralName.end()){
        auto &shortName = it->first;
        if(fullName.find(shortName) != std::string::npos){
          int type = typeGeneralName[shortName];
          if (type == static_cast<int>(ChromosomeSubstitutionModel::BASENUM)){
            chromEvolParams[i].first = static_cast<int>(allParams.getParameter(fullName).getValue());  
          }else{
            if (chromEvolParams[i].second.find(type) != chromEvolParams[i].second.end()){
              chromEvolParams[i].second[type].push_back(allParams.getParameter(fullName).getValue());

            }else{
              chromEvolParams[i].second[type];
              chromEvolParams[i].second[type].push_back(allParams.getParameter(fullName).getValue());

            }
            
          }
          break;
        }
        it ++;
      }
    }
    for (int j = 0; j < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; j++){
      uint numOfParams = LikelihoodUtils::getNumberOfParametersPerParamType(j, ChromEvolOptions::rateChangeType_);
      if (numOfParams == 0){
        if (j == static_cast<int>(ChromosomeSubstitutionModel::BASENUM)){
          chromEvolParams[i].first = IgnoreParam;

        }else{
          chromEvolParams[i].second[j];
        }
      }

    }
  }

  uint numOfRetainedModels = numOfRequiredModels-numOfModels;
  for (uint i = 1; i <= numOfRetainedModels; i++){
    chromEvolParams[i + numOfModels] = chromEvolParams[numOfModels]; 
  }
  return chromEvolParams;
}
/******************************************************************************************/

void ChromosomeTraitOptimizer::setParametersToNewTraitModel(std::map<string, double> &traitModelParams, std::shared_ptr<CharacterSubstitutionModel> traitModel, shared_ptr<IntegerFrequencySet> freqs, bool random){
  traitModel->setFrequencySet(*freqs);
  auto it = traitModelParams.begin();
  while (it != traitModelParams.end()){
    auto paramName = it->first;
    string prefix = "pi";
    if (paramName.compare(0, prefix.length(), prefix) != 0){
      if ((random) && (std::find(fixedTraitParams_.begin(), fixedTraitParams_.end(), paramName) == fixedTraitParams_.end())){
        auto value = RandomTools::giveRandomNumberBetweenTwoPoints(1e-3, 100);
        traitModel->setParameterValue(paramName, value);
      }else{ // if not random or the parameter is fixed
        traitModel->setParameterValue(paramName, traitModelParams[paramName]);
      }
    }
    it ++;
  }
}

void ChromosomeTraitOptimizer::optimizeJointLikelihood()
{
  //Go over each cycle
  // note: I start from the second point, because the first iteration has been already done.

  for (size_t i = 0; i < numOfIterations_.size(); i++){
    if (i > 0){
      cout <<"##################################" << endl;
      cout << "*********  cycle " << i << "  **************"<<endl;  

    }
    for (size_t j = 0; j < numOfPoints_[i]; j++){
      //If the number of optimization iterations is larger than zero, optimize the number of times as specified
      std::cout << "Optimizing point #" << j << "..." << std::endl;
      if (numOfIterations_[i] > 0){
        vectorOfJointLikelohoods_[j]->optimize(tolerance_, numOfIterations_[i], numOfIterationInBetween_, baseNumberUpperBound_, baseNumOptimizationMethod_, fixedTraitRootFreq_, fixedTraitParams_);
      }

    }
    sort(vectorOfJointLikelohoods_.begin(), vectorOfJointLikelohoods_.end(), compareJointLikValues);
    if (i < numOfIterations_.size()-1){
      clearVectorOfLikelihoods(vectorOfJointLikelohoods_, numOfPoints_[i+1]);
    }
    std::cout << "*** *** *** ***" << std::endl;
    std::cout << "End of cycle " << i << std::endl;
    std::cout << "Likelihoods at the end of the cycle are:" << std::endl;
    for (size_t k = 0; k < vectorOfJointLikelohoods_.size(); k++){
      std::cout << "\t " << vectorOfJointLikelohoods_[k]->getValue() << std::endl;
    } 
    std::cout << "*** *** *** ***" << std::endl;
  }
}
