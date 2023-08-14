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

void ChromosomeTraitOptimizer::setRandomTraitModel(double* mu, double* pi0){
  *mu = RandomTools::giveRandomNumberBetweenTwoPoints(1e-3, 100);
  *pi0 = RandomTools::giveRandomNumberBetweenTwoPoints(0.05, 0.95);

}


/*
The initiation of the likelihood function is either with some fixed values or with random values. If the values are random,
it makes sense to start from assigning the same model for both states.
*/

void ChromosomeTraitOptimizer::initJointLikelihood(std::vector<double> traitModelParams, std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, const PhyloTree* tree, std::map<uint, uint> baseNumberUpperBound, std::vector<double>* rootFreqsTrait, bool random, double maxParsimonyFactor){
  std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT;
  std::shared_ptr<NonHomogeneousSubstitutionProcess> nsubProT;
  std::shared_ptr<LikelihoodCalculationSingleProcess> likT;
  //bool traitWeightedFrequencies = true;
  std::shared_ptr<BranchModel> traitModel;
  std::vector<uint> modelNodesTrait;
  std::shared_ptr<FrequencySet> rootFrequenciesTrait;
  std::shared_ptr<DiscreteDistribution> rdistTrait = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));  
  Context contextT = Context();
  if (traitModel_ == "Binary"){
    const BinaryAlphabet* traitAlpha = dynamic_cast<const BinaryAlphabet*>(vscTrait_->getAlphabet());
    double mu;
    double pi0;
    if (!random){
      mu = traitModelParams[0];
      pi0 = traitModelParams[1];
    }else{
      setRandomTraitModel(&mu, &pi0);

    }
    auto twoParamModel = std::make_shared<TwoParameterBinarySubstitutionModel>(traitAlpha,mu,pi0);
    traitModel = static_pointer_cast<BranchModel>(twoParamModel);
    subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);

    
    // if (rootFreqsTrait){
      
    //   if (fixedTraitRootFreq_){
    //     std::shared_ptr<FixedFrequencySet> rootFreqsFixedTrait = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(twoParamModel->getStateMap(), false)), *rootFreqsTrait);
    //     rootFrequenciesTrait = static_pointer_cast<FrequencySet>(rootFreqsFixedTrait);

    //   }else{
    //     if (random){
    //       setRandomRootFreqs(rootFrequenciesTrait, traitModel_);
    //       // do something
    //     }else{
    //       std::shared_ptr<FullFrequencySet> estimatedRootFreqs = std::make_shared<FullFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(twoParamModel->getStateMap(), false)));
    //       rootFrequenciesTrait = static_pointer_cast<FrequencySet>(estimatedRootFreqs);

    //     }

    //   }
      
      
    //   subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree, rootFrequenciesTrait);
    //   traitWeightedFrequencies = false;        

    // }else{
    //   subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);

    // }
    getTraitNodes(modelNodesTrait);
    
    subProT->addModel(std::shared_ptr<TwoParameterBinarySubstitutionModel>(twoParamModel->clone()), modelNodesTrait);
    nsubProT = std::shared_ptr<NonHomogeneousSubstitutionProcess>(subProT->clone());

  }
  Context* context = new Context();
  likT = std::make_shared<LikelihoodCalculationSingleProcess>(contextT, *vscTrait_->clone(), *nsubProT);
  auto phyloT = std::make_shared<SingleProcessPhyloLikelihood>(contextT, likT);
  phyloT->getValue(); // calculate now because it will be needed for stochastic mapping

  // get expected tree for the chromosome model
  StochasticMapping* stm = new StochasticMapping(likT, numberOfStochasticMappings_);
  auto treeChr = createExpectedMappingHistory(numberOfStochasticMappings_, stm);
  // creating the chromosome number model
  auto rdistC = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
  std::map<uint, vector<uint>> mapModelNodesIds = getNodesForEachModel(treeChr, stm);
  delete stm;
  std::shared_ptr<ParametrizablePhyloTree> parTreeChr;
  std::vector<std::shared_ptr<ChromosomeSubstitutionModel>> models;
  SubstitutionProcess* subProcess;
  if (random){
    subProcess = LikelihoodUtils::setRandomChromosomeSubstitutionModel(treeChr.get(), vscChr_, alphabetChr_, baseNumberUpperBound_, mapModelNodesIds, modelParams, 2, maxParsimonyFactor, fixedParams_, &sharedParams_, weightedFreqs_, &models, parTreeChr);

  }else{
    subProcess = LikelihoodUtils::setChromosomeSubstitutionModel(treeChr.get(), vscChr_, alphabetChr_, baseNumberUpperBound_, mapModelNodesIds, modelParams, 2, &sharedParams_, weightedFreqs_, &models, parTreeChr);

  }
  auto parTree2 =  std::make_shared<ParametrizablePhyloTree>(*treeChr);
  

  auto modelColl=std::make_shared<SubstitutionProcessCollection>();
  
  modelColl->addModel(traitModel, 1);
  modelColl->addModel(models[0], 2);
  modelColl->addModel(models[1], 3);
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
  std::vector<uint> vP2m1 = mapModelNodesIds[1];
  std::vector<uint> vP2m2 = mapModelNodesIds[2];
  mModBr2[2]=vP2m1;
  mModBr2[3]=vP2m2;

  modelColl->addSubstitutionProcess(2, mModBr2, 2, 2);

  // Now setting the joint likelihood object
    // Likelihoods
  SubstitutionProcess* sP1c=nsubProT->clone();
  SubstitutionProcess* sP2c=subProcess->clone(); // no need to clone again- it was already cloned.
  delete subProcess;
  auto pc(std::make_shared<PhyloLikelihoodContainer>(*context, *modelColl));



  auto lik1_j = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *(vscTrait_->clone()), *sP1c);

  pc->addPhyloLikelihood(1, new SingleProcessPhyloLikelihood(*context, lik1_j));
    
  auto lik2_j = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *(vscChr_->clone()), *sP2c, weightedFreqs_);

  pc->addPhyloLikelihood(2, new SingleProcessPhyloLikelihood(*context, lik2_j));
  //   JointTraitChromosomeLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool expectedHistory, bool weightedFrequencies, size_t numOfMappings, std::string traitModel, std::vector<int> &rateChangeType, VectorSiteContainer* chromosomeVsc, std::vector<unsigned int> &baseNumberCandidates, bool inCollection = true);
  JointTraitChromosomeLikelihood* jl = new JointTraitChromosomeLikelihood(*context, pc, true, true, numberOfStochasticMappings_, traitModel_, ChromEvolOptions::rateChangeType_, vscChr_->clone(), baseNumberCandidates_);
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

void ChromosomeTraitOptimizer::initMultipleLikelihoodPoints(std::vector<double> &traitModelParams, std::map<uint, std::pair<int, std::map<int, vector<double>>>> &modelParams, const PhyloTree* tree, std::map<uint, uint> baseNumberUpperBound, std::vector<double>* rootFreqsTrait){
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
          initJointLikelihood(traitModelParams, modelParams, tree, baseNumberUpperBound, rootFreqsTrait, false, 0);
        }else{
          double factor;
          if (n > 1){
            factor = parsimonyBound_ * (1+(0.1*(double)n));

          }else{
            factor = parsimonyBound_ * static_cast<double>(n);
          }
          
          initJointLikelihood(traitModelParams, modelParams, tree, baseNumberUpperBound, rootFreqsTrait, true, factor);
    
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
  size_t numOfActualParams = numOfAllParams-fixedTraitParams_.size()-fixedChromosomeParams;
  return numOfActualParams;
}
/*********************************************************************************************/
size_t ChromosomeTraitOptimizer::getNumberOfParametersNull(){
  size_t fixedChromosomeParams = LikelihoodUtils::getNumberOfFixedParams(optimizedChromosomeLikelihood_, fixedParams_);
  auto params = vectorOfLikelihoodsTrait_[0]->getSubstitutionModelParameters();
  vector<string> parameterNames;
  LikelihoodUtils::fixSuffixForJointLikParamNames(params, parameterNames);
  std::vector<string> traitParamNames;
  fillTraitParameters(traitModel_, parameterNames, traitParamNames);
  size_t totalNumOfTraitParams = traitParamNames.size();
  size_t totalNumberOfChromosomeParameters = optimizedChromosomeLikelihood_->getSubstitutionModelParameters().size();
  size_t numOfActualParams = totalNumOfTraitParams + totalNumberOfChromosomeParameters - fixedChromosomeParams - fixedTraitParams_.size();
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
        optimizeTraitModel<SingleProcessPhyloLikelihood>(vectorOfLikelihoodsTrait_[j], tolerance_, numOfIterations_[i], parameterNames, fixedTraitParams_);
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
void ChromosomeTraitOptimizer::initTraitLikelihoods(double mu, double pi0){
    size_t index = min((int)numOfPoints_.size()-1, 1);
    vectorOfLikelihoodsTrait_.reserve(numOfPoints_[index]);
    // If base number is one of the parameters
    cout <<"##################################" << endl;
    cout << "*********  cycle 0  **************"<<endl;  
    for (size_t n = 0; n < numOfPoints_[0]; n++){
        std::cout << "Starting cycle with Point #" << n <<"...."<<endl;
        if (n == 0){
          initTraitLikelihood(mu, pi0, false);
        }else{
          initTraitLikelihood(mu, pi0, true);
        }   
    }

    sort(vectorOfLikelihoodsTrait_.begin(), vectorOfLikelihoodsTrait_.end(), LikelihoodUtils::compareLikValues);

}
/*********************************************************************************************************/
void ChromosomeTraitOptimizer::initTraitLikelihood(double mu, double pi0, bool random){
  std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree_);
  const BinaryAlphabet* traitAlpha = dynamic_cast<const BinaryAlphabet*>(vscTrait_->getAlphabet());
  std::shared_ptr<DiscreteDistribution> rdistTrait = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
  if (random){
    setRandomTraitModel(&mu, &pi0);
  } 
  auto twoParamModel = std::make_shared<TwoParameterBinarySubstitutionModel>(traitAlpha,mu,pi0);
  //traitModel = static_pointer_cast<BranchModel>(twoParamModel);
  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);
  std::vector<uint> modelNodesTrait;
  getTraitNodes(modelNodesTrait);
      
  subProT->addModel(std::shared_ptr<TwoParameterBinarySubstitutionModel>(twoParamModel->clone()), modelNodesTrait);
  //std::shared_ptr<NonHomogeneousSubstitutionProcess> nsubProT = std::shared_ptr<NonHomogeneousSubstitutionProcess>(subProT->clone());

  Context* context = new Context();
  auto likT = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vscTrait_->clone(), *(subProT->clone()));
  auto phyloT = new SingleProcessPhyloLikelihood(*context, likT);
  vectorOfLikelihoodsTrait_.push_back(phyloT);
  std::cout << "Initial likelihood is: " << phyloT->getValue() << std::endl;
  std::cout << "Model parameters are: " << std::endl;
  auto substitutionModelParams = phyloT->getSubstitutionModelParameters();
  for (size_t i = 0; i < substitutionModelParams.size(); i++){
    std::cout << "\t" << substitutionModelParams[i].getName() << " value is " << substitutionModelParams[i].getValue() << std::endl;
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
        vectorOfJointLikelohoods_[j]->optimize(tolerance_, numOfIterations_[i], numOfIterationInBetween_, baseNumberUpperBound_, baseNumOptimizationMethod_);
      }

    }
    sort(vectorOfJointLikelohoods_.begin(), vectorOfJointLikelohoods_.end(), compareJointLikValues);
    if (i < numOfIterations_.size()-1){
      clearVectorOfLikelihoods(vectorOfJointLikelohoods_, numOfPoints_[i+1]);
    }
    for (size_t k = 0; k < vectorOfJointLikelohoods_.size(); k++){
      std::cout << "Likelihoods at the end of the cycle are: " << vectorOfJointLikelohoods_[k]->getValue() << std::endl;
    }
  }
}
