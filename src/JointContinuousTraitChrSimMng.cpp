#include "JointContinuousTraitChrSimMng.h"

using namespace std;
using namespace bpp;

void JointContinuousTraitChrSimMng::init(){
    std::ifstream file(ChromEvolOptions::modelParametersCsv_);
    std::string line;
    size_t lineNum = 0;
    std::unordered_map<size_t, string> headerMap; // Map to store column index and header name
    std::unordered_map<string, vector<double>> modelParameters; // Map to store model parameters
    

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        vector<string> row;
        size_t colNum = 0;

        while (std::getline(ss, cell, ',')) {
            row.push_back(cell);
            if (lineNum == 0) {
                headerMap[colNum] = cell; // Store header name with column index
            } else {
                string headerName = headerMap[colNum];
                modelParameters[headerName].push_back(std::stod(cell)); // Convert cell to double and store in model parameters
            }
            colNum ++;
        }
        if (row.size() > 0){
            lineNum++;
        }

    }
    numOfSimulations_ = lineNum - 1; // Subtract 1 to account for header row
    updateModelParameters(modelParameters, "gain", gain_, numOfSimulations_);
    updateModelParameters(modelParameters, "loss", loss_, numOfSimulations_);
    updateModelParameters(modelParameters, "dupl", dupl_, numOfSimulations_);
    updateModelParameters(modelParameters, "demi", demi_, numOfSimulations_);
    updateModelParameters(modelParameters, "baseNumR", baseNumR_, numOfSimulations_);
    updateModelParameters(modelParameters, "baseNum", baseNum_, numOfSimulations_);
    //updateModelParameters(modelParameters, "minChrNum", minChrNum_, numOfSimulations_);
    //updateModelParameters(modelParameters, "maxChrNum", maxChrNum_, numOfSimulations_);
    //updateModelParameters(modelParameters, "maxBaseNumTransition", maxBaseNumTransition_, numOfSimulations_);
    updateModelParameters(modelParameters, "mu", mu_, numOfSimulations_);
    updateModelParameters(modelParameters, "sigma", sigma_, numOfSimulations_);
    updateModelParameters(modelParameters, "rootChromState", rootChromState_, numOfSimulations_);
    
}

std::vector<string> JointContinuousTraitChrSimMng::getMatchingParamNames(const std::unordered_map<string, vector<double>>& modelParameters, const string& paramName) {
    std::vector<std::string> matchingParamNames;
    for (const auto& specificName : modelParameters) {
        if (specificName.first.find(paramName) != std::string::npos) {
            matchingParamNames.push_back(specificName.first);
        }
    }
    return matchingParamNames;
}


void JointContinuousTraitChrSimMng::updateModelParameters(std::unordered_map<string, vector<double>>& modelParameters, const string& paramName, vector<vector<double>>& targetVector, size_t &numOfSimulations) {
    std::vector<std::string> matchingParamNames = getMatchingParamNames(modelParameters, paramName);
    targetVector.resize(numOfSimulations);
    for (size_t i = 0; i < numOfSimulations; ++i) {
        for (size_t j = 0; j < matchingParamNames.size(); ++j) {
            string specificParamName =  paramName + std::to_string(j);
            targetVector[i].push_back(modelParameters[specificParamName][i]);
        }
    }
}

void JointContinuousTraitChrSimMng::updateModelParameters(std::unordered_map<string, vector<double>>& modelParameters, const string& paramName, vector<double>& targetVector, size_t &numOfSimulations) {
    targetVector.resize(numOfSimulations);
    for (size_t i = 0; i < numOfSimulations; ++i) {
        targetVector[i] = modelParameters[paramName][i];       
    }
}

void JointContinuousTraitChrSimMng::updateModelParameters(std::unordered_map<string, vector<double>>& modelParameters, const string& paramName, vector<uint>& targetVector, size_t &numOfSimulations) {
    targetVector.resize(numOfSimulations);
    for (size_t i = 0; i < numOfSimulations; ++i) {
        targetVector[i] = static_cast<uint>(modelParameters[paramName][i]);
    }
}
std::shared_ptr<ContinuousTraitModel> JointContinuousTraitChrSimMng::createTraitModel(const string& modelName, double mu, vector<double> sigma, uint minChrNum, uint maxChrNum, ChromosomeNumberDependencyFunction::FunctionType sigmaFuncType){
    if (modelName == "BrownianMotion"){
        return std::make_shared<BrownianMotionModel>(mu, sigma, minChrNum, maxChrNum, sigmaFuncType);
    }else if (modelName == "GenomeSizeModel"){
        return std::make_shared<GenomeSizeModel>(mu, sigma, minChrNum, maxChrNum, sigmaFuncType);
    }else{
        throw Exception("JointContinuousTraitChrSimMng::createTraitModel: Unsupported continuous trait model type: " + modelName);
    }
}
void JointContinuousTraitChrSimMng::simulateEvolutionForAllModels(const string &modelName)
{
    RandomTools::setSeed(static_cast<long>(ChromEvolOptions::seed_));
    for (size_t i = 0; i < numOfSimulations_; i++){
        std::shared_ptr<ContinuousTraitModel> traitModel = createTraitModel(modelName, mu_[i], sigma_[i], minChrNum_, maxChrNum_, ChromEvolOptions::traitDependencyFunctions_["sigma"]);
        std::shared_ptr<ChromosomeAlphabet> alphabet = std::make_shared<ChromosomeAlphabet>(minChrNum_, maxChrNum_);
        std::pair<int, std::map<int, vector<double>>> mapOfParams;
        mapOfParams.first = baseNum_[i];
        mapOfParams.second[static_cast<int>(ChromosomeSubstitutionModel::GAIN)] = gain_[i];
        mapOfParams.second[static_cast<int>(ChromosomeSubstitutionModel::LOSS)] = loss_[i];
        mapOfParams.second[static_cast<int>(ChromosomeSubstitutionModel::DUPL)] = dupl_[i];
        mapOfParams.second[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)] = demi_[i];
        mapOfParams.second[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)] = baseNumR_[i];
        auto chrModel = std::make_shared<ChromosomeBMSubstitutionModel>(mu_[i], 0.01, mu_[i], ChromEvolOptions::minContinousTraitValue_, ChromEvolOptions::maxContinousTraitValue_, ChromEvolOptions::minContinousTraitValue_, ChromEvolOptions::maxContinousTraitValue_, alphabet.get(), mapOfParams.second, mapOfParams.first, maxBaseNumTransition_, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_, ChromEvolOptions::demiOnlyForEven_, true); // any sigma will work here, because it immediately gets updated according to the trait model in the simulator
        auto jointSim = std::make_shared<JointEvolutionSimulator>(tree_, chrModel, traitModel, rootChromState_[i], ChromEvolOptions::seed_ + i, ChromEvolOptions::deltaT_);
        jointSim->simulateEvolution();
        std::unordered_map<string, int> leafChromosomeStates = jointSim->getSimulatedLeafChromosomeStates();
        std::unordered_map<string, double> leafTraitStates = jointSim->getSimulatedLeafTraitStates();
        vector<double> results;
        if (jointSim->isSuccessfulSimulation()){
            for (const string& speciesName : speciesNames_){
                results.push_back(leafTraitStates[speciesName]);
            }
            for (const string& speciesName : speciesNames_){
                results.push_back(static_cast<double>(leafChromosomeStates[speciesName]));
            }

        }else{
            results = vector<double>(speciesNames_.size()*2, std::numeric_limits<double>::quiet_NaN());
        }

        simulationsResults_[i] = results;

    }

}

void JointContinuousTraitChrSimMng::writeSimulationsResultsToFile(const string& outputFilePath){
    std::ofstream outputFile(outputFilePath);
    if (!outputFile.is_open()) {
        throw Exception("JointContinuousTraitChrSimMng::writeSimulationsResultsToFile(): Could not open output file: " + outputFilePath);
    }
    for (size_t i = 0; i < speciesNames_.size() * 2; i++) {
        if (i != speciesNames_.size() * 2 - 1) {
            outputFile << speciesNames_[i % speciesNames_.size()] << ",";

        }else {
            outputFile << speciesNames_[i % speciesNames_.size()] << "\n";
        }
    }
    for (size_t i = 0; i < simulationsResults_.size(); ++i) {
        for (size_t j = 0; j < simulationsResults_[i].size(); ++j) {
            if (j != simulationsResults_[i].size() - 1) {
                outputFile << simulationsResults_[i][j] << ",";
            } else {
                outputFile << simulationsResults_[i][j] << "\n";
            }
        }
    }
    outputFile.close();
}