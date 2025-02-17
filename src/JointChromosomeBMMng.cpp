#include "JointChromosomeBMMng.h"

using namespace std;
using namespace bpp;

void JointChromosomeBMMng::setTraitData(const std::string& filePath) {
    std::ifstream inputFile(filePath);
    
    if (!inputFile.is_open()) {
        throw std::runtime_error("Error: Unable to open file.");
    }

    std::string line;
    std::string currentSpecies;
    minTraitStateInData_ = std::numeric_limits<double>::infinity();
    maxTraitStateInData_ = -std::numeric_limits<double>::infinity();
    while (std::getline(inputFile, line)) {
        // Check if the line starts with '>'
        if (!line.empty() && line[0] == '>') {
            currentSpecies = line.substr(1); // Remove the '>' character
        } else if (!line.empty()) {
            // Convert the line to a double and map it to the current species
            try {
                double traitValue = std::stod(line);
                if (traitValue < minTraitStateInData_){
                    minTraitStateInData_ = traitValue;
                }
                if (traitValue > maxTraitStateInData_){
                    maxTraitStateInData_ = traitValue;
                }
                traitData_[currentSpecies] = traitValue;
            } catch (const std::invalid_argument& e) {
                throw std::runtime_error("Error: Invalid trait value in file.");
            }
        }
    }

    inputFile.close();
}
/******************************************************************************/
std::unordered_map<string, double> JointChromosomeBMMng::getTraitData(const std::string& filePath){
    std::ifstream inputFile(filePath);
    std::unordered_map<string, double> traitData;
    if (!inputFile.is_open()) {
        throw std::runtime_error("Error: Unable to open file.");
    }

    std::string line;
    std::string currentSpecies;
    while (std::getline(inputFile, line)) {
        // Check if the line starts with '>'
        if (!line.empty() && line[0] == '>') {
            currentSpecies = line.substr(1); // Remove the '>' character
        } else if (!line.empty()) {
            // Convert the line to a double and map it to the current species
            try {
                double traitValue = std::stod(line);
                traitData[currentSpecies] = traitValue;
            } catch (const std::invalid_argument& e) {
                throw std::runtime_error("Error: Invalid trait value in file.");
            }
        }
    }
    inputFile.close();
    return traitData;
}
/******************************************************************************/
void JointChromosomeBMMng::determineMinMaxStates(){
    minTraitValue_ = minTraitStateInData_;
    maxTraitValue_ = maxTraitStateInData_;
    auto bmLikelihood = std::make_shared<BrownianMotionLikelihood>(tree_, traitData_);
    double muMLE = bmLikelihood->getMuMLE();
    double sigmaMLE = bmLikelihood->getSigmaMLE();
    vector<double> muValues;
    muValues.push_back(minTraitValue_);
    muValues.push_back(maxTraitValue_);
    for (size_t i = 0; i < muValues.size(); i++){
        std::shared_ptr<BrownianMotionAncestralReconstruction> bmAncestral = std::make_shared<BrownianMotionAncestralReconstruction>(muMLE, sigmaMLE, muValues[i], sigmaMLE, traitData_, tree_);
        bmAncestral->reconstructAncestralStates();
        auto reconstructedStates = bmAncestral->getAncestralStates();
        auto it = reconstructedStates.begin();
        while (it != reconstructedStates.end()){
            auto currState = it->second;
            if (currState > maxTraitValue_){
                maxTraitValue_ = currState;

            }else if(currState < minTraitValue_){
                minTraitValue_ = currState;
            }
            it ++;
        }

    }


}


/******************************************************************************/

std::shared_ptr<NonHomogeneousSubstitutionProcess> JointChromosomeBMMng::setHeterogeneousTraitDependentModel(std::shared_ptr<ChromosomeBMSubstitutionModel> chrModel, std::shared_ptr<ParametrizablePhyloTree> parTree, const std::unordered_map<uint, double> &statesInBranches, const PhyloTree* tree, std::shared_ptr<FrequencySet> &freqs){
    auto rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    std::shared_ptr<NonHomogeneousSubstitutionProcess> modelSet = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdist->clone()), parTree, freqs);
    // assigning for each branch a specific model (because of the branch specific state, which is the only parameter that differs between the different models)
    vector<unsigned int> nodesIds = tree->getAllNodesIndexes();
    unsigned int rootId = tree->getRootIndex();
    nodesIds.erase(std::remove(nodesIds.begin(), nodesIds.end(), rootId), nodesIds.end());
    std::sort(nodesIds.begin(), nodesIds.end());
    for (auto &nodeId : nodesIds){
        auto model = shared_ptr<ChromosomeBMSubstitutionModel>(chrModel->clone());
        model->setParameterValue("state", statesInBranches.at(nodeId));
        modelSet->addModel(model, vector<unsigned int>(1, nodeId));
    }
    return modelSet;

}

std::shared_ptr<NonHomogeneousSubstitutionProcess> JointChromosomeBMMng::initHeterogeneousModel(std::shared_ptr<ChromosomeBMSubstitutionModel> chrModel, std::shared_ptr<ParametrizablePhyloTree> parTree) const{
    auto rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    std::shared_ptr<NonHomogeneousSubstitutionProcess> modelSet = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdist->clone()), parTree);
    // assigning for each branch a specific model (because of the branch specific state, which is the only parameter that differs between the different models)
    vector<unsigned int> nodesIds = tree_->getAllNodesIndexes();
    unsigned int rootId = tree_->getRootIndex();
    nodesIds.erase(std::remove(nodesIds.begin(), nodesIds.end(), rootId), nodesIds.end());
    std::sort(nodesIds.begin(), nodesIds.end());
    for (auto &nodeId : nodesIds){
        modelSet->addModel(shared_ptr<ChromosomeBMSubstitutionModel>(chrModel->clone()), vector<unsigned int>(1, nodeId));
    }
    // Aliasing the global parameters
    auto params = chrModel->getParameters();
    auto paramNames = params.getParameterNames(); // this should be without suffices
    size_t numberOfModels = nodesIds.size();
    for (auto &paramName : paramNames){
        // alias a parameter only if it it not the "state" parameter
        if (paramName.find("state") == std::string::npos){
            for (size_t i = 2; i < numberOfModels + 1; i++){
                modelSet->aliasParameters(paramName + "_1", paramName + "_" + TextTools::toString(i));
            }
        }
    }
    return modelSet;
}

/******************************************************************************/

JointPhyloChromosomeBMLikelihood* JointChromosomeBMMng::initJointBMModel(std::shared_ptr<ParametrizablePhyloTree> parTree, std::shared_ptr<NonHomogeneousSubstitutionProcess> modelSet){
    SubstitutionProcess* subProcess = modelSet->clone();
    Context* context = new Context();
    size_t numberOfModels = modelSet->getNumberOfModels();
    auto modelColl=std::make_shared<SubstitutionProcessCollection>();
    
    for (size_t i = 0; i < numberOfModels; i++){
        auto model = std::dynamic_pointer_cast<const ChromosomeBMSubstitutionModel>(modelSet->getModel(i+1));
        modelColl->addModel(shared_ptr<ChromosomeBMSubstitutionModel>(model->clone()), i+1);

    }
    std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    
    modelColl->addDistribution(rdist, 1);

    modelColl->addTree(parTree, 1);
    map<size_t, Vuint> mModBr1;

    for (size_t i = 0; i < numberOfModels; i++){
      mModBr1[i+1]= modelSet->getNodesWithModel(i+1);
    }
    modelColl->addSubstitutionProcess(1, mModBr1, 1, 1);
    auto pc(std::make_shared<PhyloLikelihoodContainer>(*context, *modelColl));

      
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *(chromosomeData_->clone()), *subProcess, true);

    pc->addPhyloLikelihood(1, new SingleProcessPhyloLikelihood(*context, lik));
    return new JointPhyloChromosomeBMLikelihood(*context, pc, lik, tree_, traitData_);
}

/******************************************************************************/


void JointChromosomeBMMng::separateTraitAndChromosomeParamNames(JointPhyloChromosomeBMLikelihood* lik, vector<string> &traitParamNames, vector<string> &chromosomeParamNames) const{
    auto parameters = lik->getParameters();
    auto paramsNames = parameters.getParameterNames();
    string prefix = "Chromosome";
    for (size_t i = 0; i < paramsNames.size(); i++){
        if (paramsNames[i].size() >= prefix.size() && paramsNames[i].compare(0, prefix.size(), prefix) == 0){
            if (paramsNames[i].find("state") != std::string::npos){
                continue;
            }else if (paramsNames[i].find("mu") != std::string::npos){
                traitParamNames.push_back(paramsNames[i]);

            }else if (paramsNames[i].find("sigma") != std::string::npos){
                traitParamNames.push_back(paramsNames[i]);
            }else{
                chromosomeParamNames.push_back(paramsNames[i]);
            }
        }
    }
}

/******************************************************************************/

void JointChromosomeBMMng::optimizePointInCycle(double tol, uint numOfIterations, uint numberOfIterationsPerOneOptimization, JointPhyloChromosomeBMLikelihood* likPhyloObj){    
    // setting maps of parameter type and the corresponding parameters, and vice versa
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model
    LikelihoodUtils::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, chromosomeParamNames_, 0);
    double prevLikelihood;
    vector<string> fixed_params_names;
    for (size_t i = 0; i< numOfIterations; i++){
        std::cout << "***   ***  Outer Iteration #" << i << " ***  ***" <<  std::endl;
        likPhyloObj->setTraitOptimization();
        //std::cout << "Cycle " << i << std::endl;
        prevLikelihood = likPhyloObj->getValue();
        std::cout << "Joint likelihood value (before trait) " << prevLikelihood << std::endl;
        optimizeTraitModel(likPhyloObj, tol, numberOfIterationsPerOneOptimization, traitParamNames_, fixed_params_names);
        std::cout << "Joint likelihood value (after trait) " << likPhyloObj->getValue() << std::endl;
        std::cout << "*** Optimizing chromosome parameters ... " << std::endl;
        likPhyloObj->setChromosomeOptimization();
        optimizeChromosomeModel(tol, numberOfIterationsPerOneOptimization, chromosomeParamNames_, typeWithParamNames, paramNameAndType, likPhyloObj);
        std::cout << "Joint likelihood value (after chromosome) " << likPhyloObj->getValue() << std::endl;
        std::cout << "** End of Iteration **" << std::endl;
        if (std::abs(prevLikelihood - likPhyloObj->getValue()) < tol){
            break;
        }
    }   
}

/******************************************************************************/


void JointChromosomeBMMng::optimizeChromosomeModel(double tol, uint numOfIterations, std::vector<string> &chromosomeParamNames, std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>> &paramNameAndType, JointPhyloChromosomeBMLikelihood* likPhyloObj){
    std::map<uint, vector<int>> fixedParams;
    fixedParams[1];
    ChromosomeNumberOptimizer* chrOptimizer = new ChromosomeNumberOptimizer(likPhyloObj,
                dynamic_cast<const ChromosomeAlphabet*>(chromosomeData_->getAlphabet()), chromosomeData_, baseNumberUpperBound_, fixedParams);
    vector<uint> numOfPoints = {1};
    vector<uint> numOfItarationsVec = {numOfIterations};
    chrOptimizer->initOptimizer(numOfPoints, numOfItarationsVec, numOfPoints, numOfItarationsVec, ChromEvolOptions::optimizationMethod_, baseNumOptimizationMethod_,
        tol, ChromEvolOptions::standardOptimization_, ChromEvolOptions::BrentBracketing_, 
        ChromEvolOptions::probsForMixedOptimization_);
    chrOptimizer->optimizeModelParameters<JointPhyloChromosomeBMLikelihood>(likPhyloObj, tol, numOfIterations, baseNumberCandidates_, &(ChromEvolOptions::sharedParameters_), &(ChromEvolOptions::fixedParams_), 0, baseNumberUpperBound_, &chromosomeParamNames);
    delete chrOptimizer;
    

}

/******************************************************************************/

void JointChromosomeBMMng::clearVectorOfLikelihoods(vector<JointPhyloChromosomeBMLikelihood*> &vectorOLikelihoods, size_t new_size){
    while(vectorOLikelihoods.size() > new_size){
        JointPhyloChromosomeBMLikelihood* lik_to_del = vectorOLikelihoods.back(); 
        vectorOLikelihoods.pop_back();
        if (lik_to_del){
            delete lik_to_del;

        }
    }
    vectorOLikelihoods.shrink_to_fit();
}

/******************************************************************************/

double JointChromosomeBMMng::getRoughSigmaEstimator(){
    double average = 0;
    vector<double> states;
    for (const auto& it : traitData_) {
        states.push_back(it.second);
        average += it.second;
    }
    average /= static_cast<double>(states.size());
    double variance = 0;
    for (auto &state : states){
        variance += ((state - average) *  (state - average));
    }
    variance /= static_cast<double>(states.size());
    return variance;
}

/******************************************************************************/

void JointChromosomeBMMng::initJointPhyloChromosomeBMLikelihood(ChromosomeAlphabet* alphabet, pair<int, std::map<int, std::vector<double>>> &modelsParams, ChromosomeSubstitutionModel::rootFreqType freqType,
    vector<int> &rateChangeType,
    bool demiOnlyForEven,
    double sigma, double mu, double estimatedSigma, double parsimonyBound, bool random)
{
    std::shared_ptr<ChromosomeBMSubstitutionModel> chrModel;
    if (random){
        vector<int> fixedParams;
        chrModel = std::shared_ptr<ChromosomeBMSubstitutionModel>(ChromosomeBMSubstitutionModel::initBMRandomModel(alphabet, modelsParams.first, modelsParams.second, baseNumberUpperBound_[1], freqType,
         rateChangeType, fixedParams, demiOnlyForEven, parsimonyBound, minTraitValue_, maxTraitValue_, minTraitStateInData_, maxTraitStateInData_, estimatedSigma));

    }else{
        chrModel = std::make_shared<ChromosomeBMSubstitutionModel>(mu, sigma, 0, minTraitValue_, maxTraitValue_, minTraitStateInData_, maxTraitStateInData_, alphabet, modelsParams.second, modelsParams.first, baseNumberUpperBound_[1], freqType, rateChangeType, demiOnlyForEven);

    }
    auto parTree =  std::make_shared<ParametrizablePhyloTree>(*tree_);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> modelSet = initHeterogeneousModel(chrModel, parTree);
    JointPhyloChromosomeBMLikelihood* likObj = initJointBMModel(parTree, modelSet);
    likObj->getInitialLikelihood();
    vectorOfJointLikelohoods_.push_back(likObj);
    std::cout << "Initial likelihood is " << likObj->getValue() << std::endl;
    auto params = likObj->getSubstitutionModelParameters();
    std::cout <<"Parameters are:" << std::endl;
    for (size_t i = 0; i < params.size(); i++){
        auto paramName = params[i].getName();
        if (paramName.find("state") != std::string::npos){
            continue;
        }
        std::cout <<"\t" << paramName << ": " << params[i].getValue() << std::endl;

    }
}

/******************************************************************************/
void JointChromosomeBMMng::initMultipleLikelihoodPoints(ChromosomeAlphabet* alphabet, pair<int, std::map<int, std::vector<double>>> &modelsParams, ChromosomeSubstitutionModel::rootFreqType freqType,
    vector<int> &rateChangeType,
    bool demiOnlyForEven,
    double sigma, double mu, double parsimonyBound, int &seed)
{
    double estimatedSigma = getRoughSigmaEstimator();
    std::cout << "******************************************************" << std::endl;
    std::cout << "*  Optimizing joint likelihood function  *" << std::endl;
    std::cout << "******************************************************" << std::endl;
    size_t index = min((int)numOfPoints_.size()-1, 1);
    vectorOfJointLikelohoods_.reserve(numOfPoints_[index]);
    // If base number is one of the parameters
    cout <<"##################################" << endl;
    cout << "*********  cycle 0  **************"<<endl;
    int count_num_of_trials; 
    int max_trials = 10;
    RandomTools::setSeed(static_cast<long>(seed));
    for (size_t n = 0; n < numOfPoints_[0]; n++){
      count_num_of_trials = 0;
        
      std::cout << "Starting cycle with Point #" << n <<"...."<<endl;
      if (n == 0){
        initJointPhyloChromosomeBMLikelihood(alphabet, modelsParams, freqType, rateChangeType, demiOnlyForEven, sigma, mu, estimatedSigma, parsimonyBound, false);
        separateTraitAndChromosomeParamNames(vectorOfJointLikelohoods_[0], traitParamNames_, chromosomeParamNames_);

      }else{
        double factor;
        if (n > 1){
          factor = parsimonyBound * (1+(0.1*(double)n));

        }else{
          factor = parsimonyBound * static_cast<double>(n);
        }

        while (count_num_of_trials < max_trials) {
          try {
            initJointPhyloChromosomeBMLikelihood(alphabet, modelsParams, freqType, rateChangeType, demiOnlyForEven, sigma, mu, estimatedSigma, parsimonyBound, true);
            count_num_of_trials ++;
            break; // If successful, exit the loop
          }catch (const std::runtime_error& e) {
            count_num_of_trials++;
            if (count_num_of_trials >= max_trials) {
              throw std::runtime_error(std::string("Mapping failure after ") + std::to_string(count_num_of_trials) + " attempts: " + e.what());
            }
          }
        }          
      }
      if (numOfPoints_.size() > 1){
        removeUnnecessaryLikObject(vectorOfJointLikelohoods_, numOfPoints_[1]);
      }   
    }
    //sort(vectorOfJointLikelohoods_.begin(), vectorOfJointLikelohoods_.end(), compareJointLikValues);
}
/******************************************************************************/
void JointChromosomeBMMng::removeUnnecessaryLikObject(vector<JointPhyloChromosomeBMLikelihood*> &vectorOLikelihoods, size_t numberOfBestPoints){
    if (vectorOLikelihoods.size() <= numberOfBestPoints){
        return;
    }
    sort(vectorOLikelihoods.begin(), vectorOLikelihoods.end(), compareJointLikValues);
    clearVectorOfLikelihoods(vectorOLikelihoods, numberOfBestPoints);

}

/******************************************************************************/

void JointChromosomeBMMng::optimizeJointLikelihood()
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
        optimizePointInCycle(tolerance_, numOfIterations_[i], numOfIterationsPerStep_, vectorOfJointLikelohoods_[j]);
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

/******************************************************************************/
void JointChromosomeBMMng::optimizeTraitModel(JointPhyloChromosomeBMLikelihood* lik, double tol, uint numOfIterations, std::vector<string> &traitParamNames, std::vector<string> &fixedParamsTrait){
    DerivableSecondOrder* f = lik;
    BrentOneDimension* optimizer = new BrentOneDimension(f);
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);
    // setting bracketing for Brent optimization
    optimizer->setBracketing(BrentOneDimension::BRACKET_SIMPLE);
    ParameterList params;
    double currentLikelihood = lik->getValue();
    double prevLikelihood;
    std::cout << "Optimizing trait model ...." << std::endl;
    for (size_t i = 0; i < numOfIterations; i++){
        std::cout << "Iteration #" << i << std::endl;
        for (size_t j = 0; j < traitParamNames.size(); j++){
            prevLikelihood = currentLikelihood;
            params = lik->getParameters();
            const string nameOfParam = traitParamNames[j];
            std::cout << "Previous value of " << nameOfParam  << " is: "+ std::to_string(params.getParameter(nameOfParam).getValue()) << std::endl;
            if (std::count(fixedParamsTrait.begin(), fixedParamsTrait.end(), nameOfParam)){
                continue;
            }
            Parameter param = params.getParameter(nameOfParam);
            if (nameOfParam.find("mu_") != std::string::npos){
                lik->setMuOptimization();
            }
            std::cout << "Parameter name is: " << nameOfParam << std::endl;
            std::shared_ptr<IntervalConstraint> bounds = std::dynamic_pointer_cast<IntervalConstraint>(lik->getParameter(nameOfParam).getConstraint());
                 
            if ((i == 1) & (numOfIterations > 2)){
                optimizer->getStopCondition()->setTolerance(tol* 2);
            }else{
                optimizer->getStopCondition()->setTolerance(tol);
            }
            optimizer->setInitialInterval(bounds->getLowerBound()+0.00001, bounds->getUpperBound()-0.00001); 
           
            optimizer->init(params.createSubList(param.getName()));
            currentLikelihood = optimizer->optimize();
            std::cout << "parameter value after optimization "+ std::to_string(lik->getParameter(param.getName()).getValue()) << std::endl;
            std::cout << "***" << std::endl;
            if (nameOfParam.find("mu_") != std::string::npos){
                lik->unsetMuOptimization();
            }
        }
        printLikelihoodParameters(lik, 1, traitParamNames);
        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
    }
    delete optimizer;

}
/******************************************************************************/

void JointChromosomeBMMng::printLikelihoodParameters(JointPhyloChromosomeBMLikelihood* lik, unsigned int optimized, vector<string> paramsNames) const{
    if (optimized == 0){
        std::cout << "Initial likelihood is: " << lik->getValue() << std::endl;
        std::cout << "Parameters are:" << std::endl;
    }else{
        std::cout << "Optimized likelihood is: " << lik->getValue() << std::endl;
        std::cout << "Optimized parameters are:"<< std::endl;
    }   
    for (int i = 0; i < (int)(paramsNames.size()); i++){

        std::cout << paramsNames[i] << " = " << lik->getParameter(paramsNames[i]).getValue() << std::endl;

             
    }
    std::cout <<  "***" << std::endl;

}

/******************************************************************************/

size_t JointChromosomeBMMng::getNumberOfParameters() const{
    return chromosomeParamNames_.size() + traitParamNames_.size();

}
/******************************************************************************/
string JointChromosomeBMMng::getParamWithoutSuffixAndPrefix(string &paramName){
    std::string suffix = "_1";
    string paramNameCorrected;
    if (paramName.size() >= suffix.size() && paramName.compare(paramName.size() - suffix.size(), suffix.size(), suffix) == 0) {
        paramNameCorrected = paramName.substr(0, paramName.size() - suffix.size());
    }
    else{
        throw Exception("JointChromosomeBMMng::getParamWithoutSuffixAndPrefix(): ERROR!! string has a different suffix");
    }
    string prefix = "Chromosome.";
    if (paramNameCorrected.find(prefix) == 0) {
        paramNameCorrected.erase(0, prefix.length());
    }
    return paramNameCorrected;

}
/******************************************************************************/
std::shared_ptr<NonHomogeneousSubstitutionProcess> JointChromosomeBMMng::setSubstitutionProcess(ValueRef <Eigen::RowVectorXd> rootFreqs, SingleProcessPhyloLikelihood* chromosomeLik, const JointPhyloChromosomeBMLikelihood* lik, std::shared_ptr<ParametrizablePhyloTree> parTree){
    if (parTree == nullptr){
        parTree =  std::make_shared<ParametrizablePhyloTree>(*tree_);

    }
    auto rootFreqsValues =  rootFreqs->getTargetValue();
    Vdouble rootFreqsBpp;
    copyEigenToBpp(rootFreqsValues, rootFreqsBpp);
    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chromosomeLik->getSubstitutionProcess().getModel(1)->getStateMap(), false)), rootFreqsBpp);
    std::shared_ptr<FrequencySet> rootFrequencies = std::shared_ptr<FrequencySet>(rootFreqsFixed->clone());
    auto rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    std::shared_ptr<NonHomogeneousSubstitutionProcess> modelSet = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdist->clone()), parTree, rootFrequencies);
    // assigning for each branch a specific model (because of the branch specific state, which is the only parameter that differs between the different models)
    vector<unsigned int> nodesIds = tree_->getAllNodesIndexes();
    unsigned int rootId = tree_->getRootIndex();
    nodesIds.erase(std::remove(nodesIds.begin(), nodesIds.end(), rootId), nodesIds.end());
    std::sort(nodesIds.begin(), nodesIds.end());
    for (size_t i = 0; i < nodesIds.size(); i++){
        auto modelToDuplicate = std::dynamic_pointer_cast<const ChromosomeBMSubstitutionModel>(chromosomeLik->getSubstitutionProcess().getModel(i+1));
        auto duplicatedModel = std::shared_ptr<ChromosomeBMSubstitutionModel>(modelToDuplicate->clone());
        for (auto &paramName: chromosomeParamNames_){
            string modelSpecificParamName = getParamWithoutSuffixAndPrefix(paramName);
            duplicatedModel->setParameterValue(modelSpecificParamName, lik->getParameter(paramName).getValue());
        }
        for (auto &paramName : traitParamNames_){
            string modelSpecificParamName = getParamWithoutSuffixAndPrefix(paramName);
            duplicatedModel->setParameterValue(modelSpecificParamName, lik->getParameter(paramName).getValue());
        }
        duplicatedModel->setParameterValue("state", lik->getParameter("Chromosome.state_" + std::to_string(i+1)).getValue()); 
        modelSet->addModel(duplicatedModel, vector<unsigned int>(1, nodesIds[i]));
    }
    return modelSet;


}

std::shared_ptr<LikelihoodCalculationSingleProcess> JointChromosomeBMMng::setFixedJointPhyloBMLikelihoodForMLAncestralReconstruction(){
    auto parTree =  std::make_shared<ParametrizablePhyloTree>(*tree_);
    auto lik = getBestLikelihoodObject();
    auto chromosomeLik = const_cast<SingleProcessPhyloLikelihood*>(dynamic_cast<const SingleProcessPhyloLikelihood*>(lik->getAbstractPhyloLikelihood(lik->getNumbersOfPhyloLikelihoods()[0])));
    ValueRef <Eigen::RowVectorXd> rootFreqs = chromosomeLik->getLikelihoodCalculationSingleProcess()->getRootFreqs();
    auto modelSet = setSubstitutionProcess(rootFreqs, chromosomeLik, lik, parTree);
    SubstitutionProcess* subProcess = modelSet->clone();
    Context* context = new Context();
    std::shared_ptr<LikelihoodCalculationSingleProcess> ancestralLik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *chromosomeData_->clone(), *subProcess, rootFreqs);
    return ancestralLik;

}
/***************************************************************************************************************************** */
std::unordered_map<uint, double> JointChromosomeBMMng::getBranchStatesForEachNode(){
    std::unordered_map<uint, double> branchStates;
    auto chromLik = dynamic_cast<const SingleProcessPhyloLikelihood*>(vectorOfJointLikelohoods_[0]->getAbstractPhyloLikelihood(vectorOfJointLikelohoods_[0]->getNumbersOfPhyloLikelihoods()[0]));
    const NonHomogeneousSubstitutionProcess* substitutionProcess = dynamic_cast<const NonHomogeneousSubstitutionProcess*>(&chromLik->getSubstitutionProcess());
    size_t numberOfModels = substitutionProcess->getNumberOfModels();
    for (size_t i = 0; i < numberOfModels; i++){
        uint nodeId = substitutionProcess->getNodesWithModel(i+1)[0];
        string paramName = "Chromosome.state_" + std::to_string(i+1);
        auto branchValue = vectorOfJointLikelohoods_[0]->getParameter(paramName).getValue();
        branchStates[nodeId] = branchValue;

    }
    return branchStates;

}