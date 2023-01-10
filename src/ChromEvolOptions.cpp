
#include "ChromEvolOptions.h"
using namespace bpp;
using namespace std;

string ChromEvolOptions::treeFilePath_;
string ChromEvolOptions::characterFilePath_;
int ChromEvolOptions::maxChrNum_;
int ChromEvolOptions::minChrNum_;
int ChromEvolOptions::numOfModels_;
double ChromEvolOptions::branchMul_;
std::vector <unsigned int> ChromEvolOptions::OptPointsNum_;
std::vector <unsigned int> ChromEvolOptions::OptIterNum_;
std::vector <unsigned int> ChromEvolOptions::OptPointsNumNextRounds_;
std::vector <unsigned int> ChromEvolOptions::OptIterNumNextRounds_;
std::map<uint, vector<double>> ChromEvolOptions::gain_;
std::map<uint, vector<double>> ChromEvolOptions::loss_;
std::map<uint, vector<double>> ChromEvolOptions::dupl_;
std::map<uint, vector<double>> ChromEvolOptions::demiDupl_;
std::map<uint, int> ChromEvolOptions::baseNum_;
std::map<uint, vector<double>> ChromEvolOptions::baseNumR_;
double ChromEvolOptions::tolerance_;
unsigned int ChromEvolOptions::maxIterations_;
bool ChromEvolOptions::maxParsimonyBound_;
unsigned int ChromEvolOptions::maxAlpha_;
unsigned int ChromEvolOptions::minAlpha_;
int ChromEvolOptions::BrentBracketing_;
bool ChromEvolOptions::standardOptimization_;
string ChromEvolOptions::optimizationMethod_;
int ChromEvolOptions::seed_;
std::vector <double> ChromEvolOptions::probsForMixedOptimization_;
string ChromEvolOptions::rootFreqs_;
string ChromEvolOptions::fixedFrequenciesFilePath_;
std::vector<int> ChromEvolOptions::rateChangeType_;
//bool ChromEvolOptions::optimizeBaseNumber_;
string ChromEvolOptions::baseNumOptimizationMethod_;
std::map<uint, std::vector<int>> ChromEvolOptions::fixedParams_;
int ChromEvolOptions::NumOfSimulations_;
int ChromEvolOptions::jumpTypeMethod_;
bool ChromEvolOptions::simulateData_;
int ChromEvolOptions::numOfDataToSimulate_;
string ChromEvolOptions::resultsPathDir_;
std::map<uint, uint> ChromEvolOptions::maxBaseNumTransition_;
double ChromEvolOptions::treeLength_;
int ChromEvolOptions::maxNumOfTrials_;
int ChromEvolOptions::minCladeSize_;
int ChromEvolOptions::maxNumOfModels_;
std::map<uint, std::vector<uint>> ChromEvolOptions::mapOfNodeIdsPerModel_;
std::map<int, vector<std::pair<uint, int>>> ChromEvolOptions::sharedParameters_;
bool ChromEvolOptions::heterogeneousModel_;
double ChromEvolOptions::deltaAICcThreshold_;
std::map<uint, std::vector<uint>> ChromEvolOptions::mapModelNodesIds_;
string ChromEvolOptions::nodeIdsFilePath_;
std::vector<uint> ChromEvolOptions::initialModelNodes_;
std::vector<string> ChromEvolOptions::globalParams_;
bool ChromEvolOptions::parallelization_;
int ChromEvolOptions::maxChrInferred_;
bool ChromEvolOptions::backwardPhase_;
bool ChromEvolOptions::forwardPhase_;
bool ChromEvolOptions::runStochasticMapping_;
size_t ChromEvolOptions::numOfStochasticMappingTrials_;
size_t ChromEvolOptions::numOfFixingMappingIterations_; // number of mapping trials
bool ChromEvolOptions::useMaxBaseTransitonNumForOpt_;
string ChromEvolOptions::modelSelectionCriterion_;
size_t ChromEvolOptions::numOfSimulatedData_;
double ChromEvolOptions::fracAllowedFailedSimulations_;
bool ChromEvolOptions::correctBaseNumber_;
size_t ChromEvolOptions::numOfRequiredSimulatedData_;

/*************************************************************************/
void ChromEvolOptions::initAllParameters(BppApplication& ChromEvol){
    initDefaultParameters();
    initParametersFromFile(ChromEvol);

}
/*************************************************************************/
void ChromEvolOptions::initDefaultParameters(){
    maxAlpha_ = 500;
    minAlpha_ = 1;
    maxChrNum_ = -10;
    minChrNum_ = 1;
    numOfModels_ = 1;
    maxIterations_ = 5;
    tolerance_ = 0.1;
    branchMul_ = 999;
    //baseNum_ = IgnoreParam;
    maxParsimonyBound_ = false;
    standardOptimization_ = false;
    BrentBracketing_ = 2;
    optimizationMethod_ = "Brent";
    seed_ = 0;
    rootFreqs_ = "weighted";
    //optimizeBaseNumber_ = false;
    baseNumOptimizationMethod_ = "Brent";
    NumOfSimulations_ = 10000;
    jumpTypeMethod_ = 0;
    simulateData_ = false;
    numOfDataToSimulate_ = 1;
    maxBaseNumTransition_[1] = 18;
    treeLength_ = 0;
    maxNumOfTrials_ = 10;
    minCladeSize_ = 2;
    maxNumOfModels_ = 1;
    heterogeneousModel_ = false; // the default is homogeneous model
    deltaAICcThreshold_ = 2;
    parallelization_ = false;
    maxChrInferred_ = maxChrNum_;
    backwardPhase_ = true;
    forwardPhase_ = true;
    runStochasticMapping_ = false;
    numOfStochasticMappingTrials_ = 1000000;
    numOfFixingMappingIterations_ = 1000;
    useMaxBaseTransitonNumForOpt_ = false;
    modelSelectionCriterion_ = "AICc";
    numOfSimulatedData_ = 1;
    fracAllowedFailedSimulations_ = 0.01;
    correctBaseNumber_ = true;
    numOfRequiredSimulatedData_ = numOfSimulatedData_;


    

}
/*************************************************************************/
// std::string ChromEvolOptions::getParamName(int type){
//     std::string paramName;
//     switch (type)
//     {
//     case ChromosomeSubstitutionModel::GAIN:
//         paramName = "gain";
//         break;
//     case ChromosomeSubstitutionModel::LOSS:
//         paramName= "loss";
//         break;
//     case ChromosomeSubstitutionModel::DUPL:
//         paramName = "dupl";
//         break;
//     case ChromosomeSubstitutionModel::DEMIDUPL:
//         paramName = "demi";
//         break;
//     case ChromosomeSubstitutionModel::BASENUM:
//         paramName = "baseNum";
//         break;
//     case ChromosomeSubstitutionModel::BASENUMR:
//         paramName = "baseNumR";
//         break;
    
//     default:
//         throw Exception("ChromEvolOptions::getParamName(): No such parameter exists!!!"); 
//         break;
//     }
//     return paramName;

// }
/*************************************************************************/
std::vector<int> ChromEvolOptions::translateStringParamsToInt(std::vector<string> &strParams){
    std::vector <int> params;
    for (size_t i = 0; i < strParams.size(); i++){
        if (strParams[i] == "gain"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::GAIN));   
        }else if (strParams[i] == "loss"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::LOSS));
        }else if (strParams[i] == "dupl"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::DUPL));
        }else if (strParams[i] == "demiPloidyR"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL));
        }else if (strParams[i] == "baseNumR"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::BASENUMR));
        }else if (strParams[i] == "baseNum"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::BASENUM));
        }else{
            throw Exception("ChromEvolOptions::translateStringParamsToInt(): No such parameter!!!");
        }
    }
    return params;
}


/*************************************************************************/
void ChromEvolOptions::initParametersFromFile(BppApplication& ChromEvol){
    maxChrNum_ = ApplicationTools::getIntParameter("_maxChrNum", ChromEvol.getParams(), maxChrNum_, "", true, 0);
    minChrNum_ = ApplicationTools::getIntParameter("_minChrNum", ChromEvol.getParams(), minChrNum_, "", true, 0);
    numOfModels_ = ApplicationTools::getIntParameter("_numOfModels", ChromEvol.getParams(), numOfModels_, "", true, 0);
    seed_ = ApplicationTools::getIntParameter("_seed", ChromEvol.getParams(), seed_, "", true, 0);
    simulateData_ = ApplicationTools::getBooleanParameter("_simulateData", ChromEvol.getParams(), simulateData_, "", true, 0);
    if (simulateData_){
        characterFilePath_ = ApplicationTools::getAFilePath("_dataFile", ChromEvol.getParams(), false, false, "", true, "none", 1);
    }else{
        characterFilePath_ = ApplicationTools::getAFilePath("_dataFile", ChromEvol.getParams(), true, true, "", true, "none", 1);
    }
    treeFilePath_ = ApplicationTools::getAFilePath("_treeFile", ChromEvol.getParams(), true, true, "", true, "none", 1);
    branchMul_ = ApplicationTools::getDoubleParameter("_branchMul", ChromEvol.getParams(), branchMul_, "", true, 0);
    maxIterations_ = (unsigned int)ApplicationTools::getIntParameter("_maxOptimizationItarations", ChromEvol.getParams(), maxIterations_, "", true, 0);
    tolerance_ = ApplicationTools::getDoubleParameter("_tolParamOptimization", ChromEvol.getParams(), tolerance_, "", true, 0);
    setModelParameters(ChromEvol);
    maxParsimonyBound_ = ApplicationTools::getBooleanParameter("_maxParsimonyBound", ChromEvol.getParams(), maxParsimonyBound_, "", true, 0);
    parallelization_ = ApplicationTools::getBooleanParameter("_parallelization", ChromEvol.getParams(), parallelization_, "", true, 0);
    standardOptimization_ = ApplicationTools::getBooleanParameter("_standardOptimization", ChromEvol.getParams(), standardOptimization_, "", true, 0);
    BrentBracketing_ = ApplicationTools::getIntParameter("_BrentBracketing", ChromEvol.getParams(), BrentBracketing_, "", true, 0);
    optimizationMethod_ = ApplicationTools::getStringParameter("_optimizationMethod", ChromEvol.getParams(), optimizationMethod_, "", true, 0);
    string defaultValForOptPointsNum = "10,3,1";
    string defaultValForOptIterNum = "0,2,5";
    string defaultValForProbsForMixedOpt = "1,0";
    string defaultValForOptPointsNumNextRounds = defaultValForOptPointsNum;
    string defaultValForOptIterNumNextRounds = defaultValForOptIterNum;
    OptPointsNum_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizePointsNum", ChromEvol.getParams(), ',', defaultValForOptPointsNum, "", true, 0);
    OptPointsNumNextRounds_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizePointsNumNextRounds", ChromEvol.getParams(), ',', defaultValForOptPointsNumNextRounds, "", true, 0);
    globalParams_ = ApplicationTools::getVectorParameter<string>("_globalParams", ChromEvol.getParams(), ',', "", "", true, 0);
    OptIterNum_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizeIterNum", ChromEvol.getParams(), ',', defaultValForOptIterNum, "", true, 0);
    OptIterNumNextRounds_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizeIterNumNextRounds", ChromEvol.getParams(), ',', defaultValForOptIterNumNextRounds, "", true, 0);
    probsForMixedOptimization_ = ApplicationTools::getVectorParameter<double>("_probsForMixedOptimization", ChromEvol.getParams(), ',', defaultValForProbsForMixedOpt, "", true, 0);
    fixedFrequenciesFilePath_ = ApplicationTools::getAFilePath("_fixedFrequenciesFilePath", ChromEvol.getParams(), false, true, "", true, "none", 0);
    nodeIdsFilePath_ =  ApplicationTools::getAFilePath("_nodeIdsFilePath", ChromEvol.getParams(), false, true, "", true, "none", 0);
    rootFreqs_ = ApplicationTools::getStringParameter("_rootFreqs", ChromEvol.getParams(), rootFreqs_, "", true, 0);
    std::string gainFunc = ApplicationTools::getStringParameter("_gainFunc", ChromEvol.getParams(), "None", "", true, 0);
    std::string lossFunc = ApplicationTools::getStringParameter("_lossFunc", ChromEvol.getParams(), "None", "", true, 0);
    std::string duplFunc = ApplicationTools::getStringParameter("_duplFunc", ChromEvol.getParams(), "None", "", true, 0);
    std::string demiDuplFunc = ApplicationTools::getStringParameter("_demiDuplFunc", ChromEvol.getParams(), "None", "", true, 0);
    std::string baseNumRFunc = ApplicationTools::getStringParameter("_baseNumRFunc", ChromEvol.getParams(), "None", "", true, 0);
    setFunctions(gainFunc, lossFunc, duplFunc, demiDuplFunc, baseNumRFunc);
    //optimizeBaseNumber_ = ApplicationTools::getBooleanParameter("_optimizeBaseNumber", ChromEvol.getParams(), optimizeBaseNumber_, "", true, 0);
    baseNumOptimizationMethod_ = ApplicationTools::getStringParameter("_baseNumOptimizationMethod", ChromEvol.getParams(), baseNumOptimizationMethod_, "", true, 0);
    NumOfSimulations_ = ApplicationTools::getIntParameter("_NumOfSimulations", ChromEvol.getParams(), NumOfSimulations_, "", true, 0);
    jumpTypeMethod_ = ApplicationTools::getIntParameter("_jumpTypeMethod", ChromEvol.getParams(), jumpTypeMethod_, "", true, 0);
    numOfDataToSimulate_ = ApplicationTools::getIntParameter("_numOfDataToSimulate", ChromEvol.getParams(), numOfDataToSimulate_, "", true, 0);
    resultsPathDir_ = ApplicationTools::getAFilePath("_resultsPathDir", ChromEvol.getParams(), false, true, "", true, "none", 0);
    uint maxBaseNumTransition = static_cast<uint>(ApplicationTools::getIntParameter("_maxBaseNumTransition", ChromEvol.getParams(), 18, "", true, 0));

    for (uint i = 1; i <= (uint)numOfModels_; i++){
        maxBaseNumTransition_[i] = maxBaseNumTransition;
    }
    treeLength_ = ApplicationTools::getDoubleParameter("_treeLength", ChromEvol.getParams(), treeLength_, "", true, 0);
    maxNumOfTrials_ = ApplicationTools::getIntParameter("_maxNumOfTrials", ChromEvol.getParams(), maxNumOfTrials_, "", true, 0);
    minCladeSize_ = ApplicationTools::getIntParameter("_minCladeSize", ChromEvol.getParams(), minCladeSize_, "", true, 0);
    maxNumOfModels_ = ApplicationTools::getIntParameter("_maxNumOfModels", ChromEvol.getParams(), maxNumOfModels_, "", true, 0);
    heterogeneousModel_ = ApplicationTools::getBooleanParameter("_heterogeneousModel", ChromEvol.getParams(), heterogeneousModel_, "", true, 0);
    deltaAICcThreshold_ = ApplicationTools::getDoubleParameter("_deltaAICcThreshold", ChromEvol.getParams(), deltaAICcThreshold_, "", true, 0);
    maxChrInferred_ = ApplicationTools::getIntParameter("_maxChrInferred", ChromEvol.getParams(), maxChrInferred_, "", true, 0);
    forwardPhase_ = ApplicationTools::getBooleanParameter("_forwardPhase", ChromEvol.getParams(), forwardPhase_, "", true, 0);
    backwardPhase_ = ApplicationTools::getBooleanParameter("_backwardPhase", ChromEvol.getParams(), backwardPhase_, "", true, 0);
    runStochasticMapping_ = ApplicationTools::getBooleanParameter("_runStochasticMapping", ChromEvol.getParams(), runStochasticMapping_, "", true, 0);
    numOfStochasticMappingTrials_ = static_cast<size_t>(ApplicationTools::getIntParameter("_numOfStochasticMappingTrials", ChromEvol.getParams(), (int)numOfStochasticMappingTrials_, "", true, 0));
    numOfFixingMappingIterations_ = static_cast<size_t>(ApplicationTools::getIntParameter("_numOfFixingMappingIterations", ChromEvol.getParams(), (int)numOfFixingMappingIterations_, "", true, 0));
    useMaxBaseTransitonNumForOpt_ = ApplicationTools::getBooleanParameter("_useMaxBaseTransitonNumForOpt", ChromEvol.getParams(), useMaxBaseTransitonNumForOpt_, "", true, 0);
    modelSelectionCriterion_ = ApplicationTools::getStringParameter("_modelSelectionCriterion", ChromEvol.getParams(), modelSelectionCriterion_, "", true, 0);
    numOfSimulatedData_ = static_cast<size_t>(ApplicationTools::getIntParameter("_numOfSimulatedData", ChromEvol.getParams(), (int)numOfSimulatedData_, "", true, 0));
    fracAllowedFailedSimulations_ = ApplicationTools::getDoubleParameter("_fracAllowedFailedSimulations", ChromEvol.getParams(), fracAllowedFailedSimulations_, "", true, 0);
    correctBaseNumber_ = ApplicationTools::getBooleanParameter("_correctBaseNumber", ChromEvol.getParams(), correctBaseNumber_, "", true, 0);
    numOfRequiredSimulatedData_ = static_cast<size_t>(ApplicationTools::getIntParameter("_numOfRequiredSimulatedData", ChromEvol.getParams(), (int)numOfRequiredSimulatedData_, "", true, 0));

}
/************************************************************************/
void ChromEvolOptions::setModelParameters(BppApplication& ChromEvol){
    std::map<uint, std::map<int, std::pair<int, vector<double>>>>  mapModelTypeValues;
    std::map<uint, std::pair<int, int>> mapModelBaseNumTypeAndVal;
    std::map<int, size_t> paramNums;
    vector<ChromosomeSubstitutionModel::paramType> modelTypeParams = {ChromosomeSubstitutionModel::GAIN, ChromosomeSubstitutionModel::LOSS, ChromosomeSubstitutionModel::DUPL, ChromosomeSubstitutionModel::DEMIDUPL, ChromosomeSubstitutionModel::BASENUM, ChromosomeSubstitutionModel::BASENUMR};
    vector<string> modelStringParams = {"_gain", "_loss", "_dupl", "_demiPloidyR", "_baseNum", "_baseNumR"};
    for(uint i = 1; i <= static_cast<uint>(numOfModels_); i++){
        for (size_t j = 0; j < modelStringParams.size(); j++){
            string paramName = modelStringParams[j] + "_"+ std::to_string(i);
            vector<string> paramNumAndValues = ApplicationTools::getVectorParameter<string>(paramName, ChromEvol.getParams(), ';', "", "", true, 0);
            vector<string> paramValues;
            int paramCat;
            if (paramNumAndValues.size() == 2){
                paramCat = std::stoi(paramNumAndValues[0]);
                if (paramNums.find(paramCat) == paramNums.end()){
                    paramNums[paramCat] = 0;
                }
                paramNums[paramCat] ++;
                StringTokenizer stoken = StringTokenizer(paramNumAndValues[1], ",");
                while (stoken.hasMoreToken()){
                    paramValues.push_back(stoken.nextToken());
                }
                if (modelTypeParams[j] == ChromosomeSubstitutionModel::BASENUM){
                    baseNum_[i] = std::stoi(paramValues[0]);
                    mapModelBaseNumTypeAndVal[i] = std::pair<int, int>();
                    mapModelBaseNumTypeAndVal[i].first = paramCat;
                    mapModelBaseNumTypeAndVal[i].second = baseNum_[i];
                    continue;
                }
            }else{
                paramCat = IgnoreParam;
                if (modelTypeParams[j] == ChromosomeSubstitutionModel::BASENUM){
                    baseNum_[i] = IgnoreParam;
                    mapModelBaseNumTypeAndVal[i] = std::pair<int, int>();
                    mapModelBaseNumTypeAndVal[i].first = paramCat;
                    mapModelBaseNumTypeAndVal[i].second = baseNum_[i];
                    continue;
                }

            }
            mapModelTypeValues[i][modelTypeParams[j]] = std::pair<int, std::vector<double>>();
            mapModelTypeValues[i][modelTypeParams[j]].first = paramCat;
            for (size_t k = 0; k < paramValues.size(); k++){
                mapModelTypeValues[i][modelTypeParams[j]].second.push_back(std::stod(paramValues[k]));
            }
            updateModelParameter(i, modelTypeParams[j], mapModelTypeValues[i][modelTypeParams[j]].second);
        }

    }
    //setSharedParametersInterModels();
    setSharedParametersPerModel(mapModelTypeValues, mapModelBaseNumTypeAndVal, paramNums);
    setFixedParameters(ChromEvol);
}
/************************************************************************/
// void ChromEvolOptions::setSharedParametersInterModels(){

//     std::vector<int> paramTypes = translateStringParamsToInt(globalParams_);
//     for (size_t i = 0; i < paramTypes.size(); i++){
//         for (size_t j = 1; j <= (size_t)numOfModels_; j++){
//             sharedParameters_[uint(i)].push_back(paramTypes[j]);
//         }
//     }
// }
/************************************************************************/
void ChromEvolOptions::updateModelParameter(uint model, int type, vector<double> paramValues){
    switch (type)
    {
    case ChromosomeSubstitutionModel::GAIN:
        gain_[model] = paramValues;
        break;
    case ChromosomeSubstitutionModel::LOSS:
        loss_[model] = paramValues;
        break;
    case ChromosomeSubstitutionModel::BASENUMR:
        baseNumR_[model] = paramValues;
        break;
    case ChromosomeSubstitutionModel::DUPL:
        dupl_[model] = paramValues;
        break;
    case ChromosomeSubstitutionModel::DEMIDUPL:
        demiDupl_[model] = paramValues;
        break;
    
    default:
        break;
    }


}
/************************************************************************/
void ChromEvolOptions::setFixedParameters(BppApplication& ChromEvol){
    for(uint i = 1; i <= static_cast<uint>(numOfModels_); i++){
        string paramName = "_fixedParams_"+ std::to_string(i);
        std::vector<string> fixedParamsStr = ApplicationTools::getVectorParameter<string>(paramName, ChromEvol.getParams(), ',', "", "", true, 0);
        fixedParams_[i] = translateStringParamsToInt(fixedParamsStr);

    }

}
/************************************************************************/
void ChromEvolOptions::setSharedParametersPerModel(std::map<uint, std::map<int, std::pair<int, vector<double>>>> mapModelTypeValues, std::map<uint, std::pair<int, int>> mapModelBaseNumTypeAndVal, std::map<int, size_t> paramNumFreqs){
    // for each parameter number I need to hold the corresponding parameter type and vector of models
    auto it = mapModelTypeValues.begin();
    while(it != mapModelTypeValues.end()){
        uint model = it->first;
        auto paramsPerType = it->second;
        auto itType = paramsPerType.begin();
        while(itType != paramsPerType.end()){
            int paramType = itType->first;
            int paramNum = paramsPerType[paramType].first;
            if ((paramNumFreqs[paramNum] < 2) || (paramNum == IgnoreParam)){
                itType++;
                continue;
            }
            std::pair<uint, int> modelAndType;
            modelAndType.first = model;
            modelAndType.second = paramType;
            sharedParameters_[paramNum].push_back(modelAndType);

            itType++;
        }
        it++;
    }  
    auto baseNumIt = mapModelBaseNumTypeAndVal.begin();
    while (baseNumIt != mapModelBaseNumTypeAndVal.end()){
         uint model = baseNumIt->first;
         int paramNum = mapModelBaseNumTypeAndVal[model].first;
         if ((paramNum == IgnoreParam) || (paramNumFreqs[paramNum] < 2)){
             baseNumIt ++;
             continue;
         }
         std::pair<uint, int> modelAndType;
         modelAndType.first = model;
         modelAndType.second = ChromosomeSubstitutionModel::BASENUM;
         sharedParameters_[paramNum].push_back(modelAndType);
         baseNumIt ++;
     }

    // auto paramNumIt = mapOfParamNumTypeModels.begin();
    // while (paramNumIt != mapOfParamNumTypeModels.end()){
    //     auto typeIterator = mapOfParamNumTypeModels[paramNumIt->first].begin();
    //     while (typeIterator != mapOfParamNumTypeModels[paramNumIt->first].end()){
    //         if (mapOfParamNumTypeModels[paramNumIt->first][typeIterator->first].size() > 1){
    //             for (size_t k = 0; k < mapOfParamNumTypeModels[paramNumIt->first][typeIterator->first].size(); k++){
    //                 sharedParameters_[typeIterator->first].push_back(mapOfParamNumTypeModels[paramNumIt->first][typeIterator->first][k]);
    //             }
    //         }
    //         typeIterator ++;
    //     }
    //     paramNumIt ++;
    // }

}

/************************************************************************/
void ChromEvolOptions::setFunctions(std::string gainFunc, std::string lossFunc, std::string duplFunc, std::string demiDuplFunc, std::string baseNumRFunc){
    for (size_t i = 0; i < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; i++){
        switch (i)
        {
        case ChromosomeSubstitutionModel::BASENUM:
            break;
        case ChromosomeSubstitutionModel::GAIN:
            rateChangeType_.push_back(getFunctionFromString(gainFunc));
            break;
        case ChromosomeSubstitutionModel::LOSS:
            rateChangeType_.push_back(getFunctionFromString(lossFunc));
            break;
        case ChromosomeSubstitutionModel::DUPL:
            rateChangeType_.push_back(getFunctionFromString(duplFunc));
            break;
        case ChromosomeSubstitutionModel::DEMIDUPL:
            rateChangeType_.push_back(getFunctionFromString(demiDuplFunc));
            break;
        case ChromosomeSubstitutionModel::BASENUMR:
            rateChangeType_.push_back(getFunctionFromString(baseNumRFunc));
            break;
   
        default:
            throw Exception("ChromEvolOptions::setFunctions: parameter not found !!!");
        }
    }

}


/*************************************************************************/
int ChromEvolOptions::getFunctionFromString(string funcStr){
    int func;
    if (funcStr == "CONST"){
        func = static_cast<int>(ChromosomeNumberDependencyFunction::CONSTANT);
        
    }else if (funcStr == "LINEAR"){
        func = static_cast<int>(ChromosomeNumberDependencyFunction::LINEAR);
    }else if (funcStr == "LINEAR_BD"){
        func = static_cast<int>(ChromosomeNumberDependencyFunction::LINEAR_BD);
    }else if (funcStr == "EXP"){
        func = static_cast<int> (ChromosomeNumberDependencyFunction::EXP);
    }else if (funcStr == "POLYNOMIAL"){
        func = static_cast<int> (ChromosomeNumberDependencyFunction::POLYNOMIAL);
    }else if (funcStr == "LOGNORMAL"){
        func = static_cast<int> (ChromosomeNumberDependencyFunction::LOGNORMAL);
    }else if (funcStr == "REVERSE_SIGMOID"){
        func = static_cast<int> (ChromosomeNumberDependencyFunction::REVERSE_SIGMOID);
    }else if (funcStr == "IGNORE"){ 
        func = static_cast<int> (ChromosomeNumberDependencyFunction::IGNORE);
    }else{
        throw Exception("ChromEvolOptions::getFunctionFromString(): No such function exists!!!");
    }
    return func;
}
/*************************************************************************/
void ChromEvolOptions::getInitialValuesForComplexParams(std::map<uint, std::pair<int, std::map<int, vector<double>>>> &mapOfParams){
    for (uint i = 1; i <= static_cast<uint>(numOfModels_); i++){
        mapOfParams[i] = std::pair<int, std::map<int, std::vector<double>>>();
        mapOfParams[i].first = baseNum_[i];
        mapOfParams[i].second[static_cast<int>(ChromosomeSubstitutionModel::GAIN)] = gain_[i];
        mapOfParams[i].second[static_cast<int>(ChromosomeSubstitutionModel::LOSS)] = loss_[i];
        mapOfParams[i].second[static_cast<int>(ChromosomeSubstitutionModel::DUPL)] = dupl_[i];
        mapOfParams[i].second[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)] = demiDupl_[i];
        mapOfParams[i].second[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)] = baseNumR_[i];

    }
   
}
/*************************************************************************/
// void ChromEvolOptions::initVectorOfChrNumParameters(vector<double>& paramVector){
//     for (size_t i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
//         switch(i){
//             case ChromosomeSubstitutionModel::BASENUM:
//                 paramVector.push_back(baseNum_);
//                 break;
//             case ChromosomeSubstitutionModel::BASENUMR:
//                 paramVector.push_back(baseNumR_);
//                 break;
//             case ChromosomeSubstitutionModel::DUPL:
//                 paramVector.push_back(constDupl_);
//                 break;
//             case ChromosomeSubstitutionModel::LOSS:
//                 paramVector.push_back(constLoss_);
//                 break;
//             case ChromosomeSubstitutionModel::GAIN:
//                 paramVector.push_back(constGain_);
//                 break;
//             case ChromosomeSubstitutionModel::DEMIDUPL:
//                 paramVector.push_back(constDemiDupl_);
//                 break;
//             case ChromosomeSubstitutionModel::LOSSR:
//                 paramVector.push_back(lossR_);
//                 break;
//             case ChromosomeSubstitutionModel::GAINR:
//                 paramVector.push_back(gainR_);
//                 break;
//             case ChromosomeSubstitutionModel::DUPLR:
//                 paramVector.push_back(duplR_);
//                 break;
//             default:
//                 throw Exception("ChromEvolOptions::initVectorOfChrNumParameters(): Invalid rate type!");
//                 break;
//         }

//     }   

// }
