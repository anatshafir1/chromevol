#include "LikelihoodUtils.h"

using namespace bpp;
using namespace std;


void LikelihoodUtils::setNodeIdsForAllModels(PhyloTree* phyltree, std::map<uint, std::vector<uint>> &mapModelNodesIds, string &path, std::vector<uint> &initialModelNodes){
    std::map<uint, std::vector<uint>> mapModelNodesTmp;
    if (path == "none"){
        auto nodes = phyltree->getAllNodes();
        for (size_t i = 0; i < nodes.size(); i++){
            uint nodeId = phyltree->getNodeIndex(nodes[i]);
            if (nodeId == phyltree->getRootIndex()){
                continue;
            }else{
                mapModelNodesIds[1].push_back(nodeId);
            }
        }
        return;
    }
    ifstream stream;
    stream.open(path.c_str());
    vector <string> lines = FileTools::putStreamIntoVectorOfStrings(stream);
    stream.close();
    PhyloTree* tree = phyltree->clone();
    std::map<uint, std::pair<uint, std::vector<uint>>> mapOfModelMRCAAndNodes;
    std::map<uint, uint> mapNodeModel;
    std::map<uint, uint> mapOriginalToAssignedModel;
    for (size_t i = 0; i < lines.size(); i ++){
        if (lines[i] == ""){
            continue;
        }
        getNodeIdsPerModelFromLine(lines[i], tree, mapOfModelMRCAAndNodes, mapOriginalToAssignedModel, initialModelNodes);

    }
    auto it_ModelNodes = mapOfModelMRCAAndNodes.begin();
    while (it_ModelNodes != mapOfModelMRCAAndNodes.end()){
        uint model = it_ModelNodes->first;
        uint nodeId = mapOfModelMRCAAndNodes[model].first;
        mapNodeModel[nodeId] = model;
        mapModelNodesTmp[model] = mapOfModelMRCAAndNodes[model].second;
        it_ModelNodes ++;
    }
    auto it = mapOfModelMRCAAndNodes.begin();
    while (it != mapOfModelMRCAAndNodes.end()){
        uint model = it->first;
        vector<uint> nodeIds = mapOfModelMRCAAndNodes[model].second;
        for (size_t i = 0; i < nodeIds.size(); i++){
            auto itNodeModel = mapNodeModel.find(nodeIds[i]);
            if (itNodeModel != mapNodeModel.end()){
                uint modelOfDescendant = mapNodeModel[nodeIds[i]];
                if (modelOfDescendant != model){
                    vector<uint> subtree = mapOfModelMRCAAndNodes[modelOfDescendant].second;
                    for (size_t j = 0; j < subtree.size(); j++){
                        auto nodeToDelIt = std::find(mapModelNodesTmp[model].begin(), mapModelNodesTmp[model].end(), subtree[j]);
                        if (nodeToDelIt != mapModelNodesTmp[model].end()){
                            mapModelNodesTmp[model].erase(std::remove(mapModelNodesTmp[model].begin(), mapModelNodesTmp[model].end(), subtree[j]),mapModelNodesTmp[model].end());

                        }
                    }
                }
            }
        }
        it ++;
    }
    
    delete tree;
    auto it_tmp = mapModelNodesTmp.begin();
    while(it_tmp != mapModelNodesTmp.end()){
        if (mapModelNodesIds.find(mapOriginalToAssignedModel[it_tmp->first]) != mapModelNodesIds.end()){
            mapModelNodesIds[mapOriginalToAssignedModel[it_tmp->first]].insert(mapModelNodesIds[mapOriginalToAssignedModel[it_tmp->first]].end(), mapModelNodesTmp[it_tmp->first].begin(), mapModelNodesTmp[it_tmp->first].end());

        }else{
            mapModelNodesIds[mapOriginalToAssignedModel[it_tmp->first]] = mapModelNodesTmp[it_tmp->first];

        }
        
        it_tmp ++;
    }

}
/*************************************************************************/
void LikelihoodUtils::getNodeIdsPerModelFromLine(string &content, PhyloTree* tree, std::map<uint, std::pair<uint, std::vector<uint>>> &modelAndNodeIds, std::map<uint,uint> &mapOriginalToAssignedModel, std::vector<uint> &initialModelNodes){
    vector<string> paramValues;
    std::regex modelPattern ("([\\d]+)");
    std::regex treePattern ("\\(([\\S]+)\\)");
    StringTokenizer stoken = StringTokenizer(content, "=");
    while (stoken.hasMoreToken()){
        paramValues.push_back(stoken.nextToken());
    }
    shared_ptr<PhyloNode> mrca_node;
    uint model;
    vector<uint> nodes;
    for (size_t i = 0; i < paramValues.size(); i++){
        std::smatch sm;
        if (i == 0){
            std::regex_search(paramValues[i], sm, modelPattern);
            model = std::stoi(sm[0]);
            if (modelAndNodeIds.find(model) != modelAndNodeIds.end()){
                auto maxModel = modelAndNodeIds.rbegin()->first;
                auto originalModel = model;
                model = maxModel+1;
                mapOriginalToAssignedModel[model] = originalModel;
            }else{
                mapOriginalToAssignedModel[model] = model;
            }

        }else{
            std::regex_search(paramValues[i], sm, treePattern);
            string speciesNonSepWithBrackets = sm[0];
            string speciesNonSep = speciesNonSepWithBrackets.substr(1, speciesNonSepWithBrackets.length()-2);

            vector<string> speciesNames;
            StringTokenizer speciesToken = StringTokenizer(speciesNonSep, ",");
            while (speciesToken.hasMoreToken()){
                speciesNames.push_back(speciesToken.nextToken());
            }
            std::map<std::string, shared_ptr<PhyloNode>> subtreeLeavesAsNodes;
            vector<shared_ptr<PhyloNode>> allLeaves = tree->getAllLeaves();
            for (size_t k = 0; k < allLeaves.size(); k++){
                subtreeLeavesAsNodes[allLeaves[k]->getName()] = allLeaves[k];
            }
            vector<shared_ptr<PhyloNode>> leaveNodesForMrca;
            for (size_t k = 0; k < speciesNames.size(); k++){
                leaveNodesForMrca.push_back(subtreeLeavesAsNodes[speciesNames[k]]);
            }
            
            mrca_node = TreeUtils::getMRCA(tree, leaveNodesForMrca);
 
            auto mrca_id = tree->getNodeIndex(mrca_node);
            initialModelNodes.push_back(mrca_id);
            auto subtreeNodes = tree->getSubtreeNodes(mrca_node);
            auto allNodeIds = tree->getNodeIndexes(subtreeNodes);
            for (size_t j = 0; j  < allNodeIds.size(); j++){
                if (allNodeIds[j] == tree->getRootIndex()){
                    continue;
                }
                nodes.push_back(allNodeIds[j]);
            }
            auto leavesUnderNode = tree->getLeavesUnderNode(mrca_node);
            std:: cout << "Model #" << mapOriginalToAssignedModel[model] << std::endl;
            for (size_t j= 0; j < leavesUnderNode.size(); j++){
                std::cout << leavesUnderNode[j]->getName() << std::endl;
            }
        }    

    }
    modelAndNodeIds[model].first = tree->getNodeIndex(mrca_node);
    modelAndNodeIds[model].second = nodes;
}
/**********************************************************************/
std::map<uint, std::vector<uint>> LikelihoodUtils::findMRCAForEachModelNodes(PhyloTree* tree, std::map<uint, vector<uint>> mapOfModelsAndNodes){
    std::map <uint, std::vector<uint>> modelWithRepresentitives;
    auto it = mapOfModelsAndNodes.begin();
    while (it != mapOfModelsAndNodes.end()){
        auto nodes = mapOfModelsAndNodes[it->first];
        for (size_t i = 0; i < nodes.size(); i++){
            uint nodeId = nodes[i];
            if (nodeId == tree->getRootIndex()){
                modelWithRepresentitives[it->first].push_back(nodeId);
                break;
            }
            auto edgeIndex =  tree->getIncomingEdges(nodeId)[0]; 
            auto fatherIndex = tree->getFatherOfEdge(edgeIndex);
            if (std::find(nodes.begin(), nodes.end(), fatherIndex) == nodes.end()){
                // if the node has no father in the list, this is the representitive in the current model
                modelWithRepresentitives[it->first].push_back(nodeId);
            }

        }
        it ++;
    }
    return modelWithRepresentitives;

}
/***************************************************************************************************************************/
std::string LikelihoodUtils::getFunctionName(int func){
    std::string functionName;
    switch(func)
    {
    case ChromosomeNumberDependencyFunction::CONSTANT:
        functionName = "CONST";
        break;
    case ChromosomeNumberDependencyFunction::LINEAR:
        functionName = "LINEAR";
        break;
    case ChromosomeNumberDependencyFunction::LINEAR_BD:
        functionName = "LINEAR_BD";
        break;
    case ChromosomeNumberDependencyFunction::EXP:
        functionName = "EXP";
        break;
    case ChromosomeNumberDependencyFunction::POLYNOMIAL:
        functionName = "POLYNOMIAL";
        break;
    case ChromosomeNumberDependencyFunction::LOGNORMAL:
        functionName = "LOGNORMAL";
        break;
    case ChromosomeNumberDependencyFunction::REVERSE_SIGMOID:
        functionName = "REVERSE_SIGMOID";
        break;
    case ChromosomeNumberDependencyFunction::IGNORE:
        functionName = "IGNORE";
        break;
    default:
        throw Exception("LikelihoodUtils::getFunctionName: parameter not found !!!");
    }
    return functionName;
}

/***************************************************************************************************************************/

string LikelihoodUtils::getStringParamName(int type){
    string strName;
    if (type == ChromosomeSubstitutionModel::BASENUM){
        strName = "baseNum";
    }else if (type == ChromosomeSubstitutionModel::BASENUMR){
        strName = "baseNumR";
    }else if (type == ChromosomeSubstitutionModel::DUPL){
        strName = "dupl";
    }else if(type == ChromosomeSubstitutionModel::DEMIDUPL){
        strName = "demi";
    }else if (type == ChromosomeSubstitutionModel::GAIN){
        strName = "gain";
    }else if (type == ChromosomeSubstitutionModel::LOSS){
        strName = "loss";
    }else{
        throw Exception("LikelihoodUtils::getStringParamName(): No such parameter exists!");
    }
    return strName;
}
// /****************************************************************************/
vector <double> LikelihoodUtils::setFixedRootFrequencies(const std::string &path, std::shared_ptr<ChromosomeSubstitutionModel> chrModel){
    ifstream stream;
    stream.open(path.c_str());
    vector <double> freqs;
    vector <string> lines = FileTools::putStreamIntoVectorOfStrings(stream);
    stream.close();
    for (size_t i = 0; i < lines.size(); i++){
        string freq_i_str = TextTools::removeSurroundingWhiteSpaces(lines[i]);
        if (freq_i_str == ""){
            continue;
        }
        double freq_i = TextTools::toDouble(freq_i_str);
        if (static_cast<unsigned int>(freqs.size()) >= chrModel->getNumberOfStates()){
            if (freq_i > 0){
                throw Exception("Invalid fixed frequencies file!");
            }

        }else{
            freqs.push_back(freq_i);
        }
        
    }
    size_t nbStates = chrModel->getNumberOfStates();
    if (freqs.size() < nbStates){
        for (size_t s = freqs.size(); s < nbStates; s++){
            freqs.push_back(0);
        }
        
    }
    if (nbStates != freqs.size()){
        throw Exception("Invalid fixed frequencies file!");
    }
    return freqs;
}
/**************************************************************************************/
void LikelihoodUtils::updateWithTypeAndCorrespondingName(std::map<std::string, int> &typeGeneralName){
    typeGeneralName["gain"] = static_cast<int>(ChromosomeSubstitutionModel::GAIN);
    typeGeneralName["loss"] = static_cast<int>(ChromosomeSubstitutionModel::LOSS);
    typeGeneralName["dupl"] = static_cast<int>(ChromosomeSubstitutionModel::DUPL);
    typeGeneralName["demi"] = static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL);
    typeGeneralName["baseNumR"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUMR);
    typeGeneralName["baseNum_"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUM);
    
}
/**********************************************************************************/
int LikelihoodUtils::getTypeOfParamFromParamName(string name){
    int type;
    std::map<std::string, int> typeGeneralName;
    LikelihoodUtils::updateWithTypeAndCorrespondingName(typeGeneralName);
    auto itParamType = typeGeneralName.begin();
    while(itParamType != typeGeneralName.end()){
        string pattern = itParamType->first;
        if (name.find(pattern) != string::npos){
            type = typeGeneralName[pattern];
            break;

        }
        itParamType ++;
    }
    return type;
}
/*******************************************************************************/
size_t LikelihoodUtils::getNumberOfFixedParams(SingleProcessPhyloLikelihood* lik, std::map<uint, vector<int>> &fixedParams){
    auto substitutionModelParams = lik->getSubstitutionModelParameters();
    auto paramNames = substitutionModelParams.getParameterNames();
    size_t numOfFixedParams = 0;
    for (size_t i = 0; i < paramNames.size(); i++){
        auto paramName = paramNames[i];
        uint model = LikelihoodUtils::getModelFromParamName(paramName);
        int type = LikelihoodUtils::getTypeOfParamFromParamName(paramName);
        auto itMap = fixedParams.find(model);
        if (itMap != fixedParams.end()){
            auto fixedTypes = fixedParams[itMap->first];
            for (size_t j = 0; j < fixedTypes.size(); j++){
                if (fixedTypes[j] == type){
                    numOfFixedParams ++;
                    break;
                }
            }
        }

    }
    return numOfFixedParams;

}
/**********************************************************************************/
uint LikelihoodUtils::getModelFromParamName(string name){
    std::regex modelPattern ("_([\\d]+)");
    std::smatch sm;
    std::regex_search(name, sm, modelPattern);
    std::string modelSuffix = sm[sm.size()-1];
    uint modelId = static_cast<uint>(stoi(modelSuffix));
    return modelId;

}
/*************************************************************************************/
uint LikelihoodUtils::getNumberOfParametersPerParamType(int paramType, vector<int> &funcTypes){
    uint numOfParams;
    if (paramType == ChromosomeSubstitutionModel::BASENUM){
        if (ChromEvolOptions::baseNum_[1] == IgnoreParam){
            numOfParams = 0;
        }else{
            numOfParams = 1;
        }
        
    }else{
        size_t startForComposite = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();
        auto funcType = funcTypes[paramType-startForComposite];
        if (funcType == ChromosomeNumberDependencyFunction::IGNORE){
            numOfParams = 0;
        }else{
            ChromosomeNumberDependencyFunction* functionOp = compositeParameter::setDependencyFunction(static_cast<ChromosomeNumberDependencyFunction::FunctionType>(funcType));
            //functionOp->setDomainsIfNeeded(alphabet_->getMin(), alphabet_->getMax());
            numOfParams = static_cast<uint>(functionOp->getNumOfParameters());
            delete functionOp;

        }

    }
    return numOfParams;


}
// /*******************************************************************************/
void LikelihoodUtils::updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>>* paramNameAndType, std::vector<std::string> namesAllParams, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::string suffix){
    std::map<string, vector<std::pair<uint, int>>> sharedParamsNames;
    if (sharedParams){
        LikelihoodUtils::createMapOfSharedParameterNames(*sharedParams, sharedParamsNames);
    }
    for (size_t i = 0; i < namesAllParams.size(); i++){
        string cleanParamName;
        if (suffix != ""){
            cleanParamName = namesAllParams[i].substr(0, namesAllParams[i].length() - suffix.length());
        }else{
            cleanParamName = namesAllParams[i];
        }
        uint modelId = LikelihoodUtils::getModelFromParamName(cleanParamName);
        int type = LikelihoodUtils::getTypeOfParamFromParamName(cleanParamName);
        //should get the type

        
        typeWithParamNames[type][modelId].push_back(namesAllParams[i]);
        if (paramNameAndType){
            (*paramNameAndType)[namesAllParams[i]] = std::pair<int, uint>(type, modelId);
        }
        if (sharedParams){
            if (sharedParamsNames.find(namesAllParams[i]) != sharedParamsNames.end()){
                for (size_t k = 0; k < sharedParamsNames[namesAllParams[i]].size(); k++ ){
                    uint model = sharedParamsNames[namesAllParams[i]][k].first;
                    int typeOfSharedParam = sharedParamsNames[namesAllParams[i]][k].second;
                    typeWithParamNames[typeOfSharedParam][model].push_back(namesAllParams[i]);

                }

            }

        }

    }

}
void LikelihoodUtils::updateSharedParameters(std::map<int, vector<std::pair<uint, int>>> &sharedParams, uint prevShift, uint numOfShifts){
    auto paramNumIt = sharedParams.begin();
    // Should Change this function !!!! TODO
    int maxNum = 0;
    vector<int> globalsAccounted;
    std::vector<int> globalParamTypes = ChromEvolOptions::translateStringParamsToInt(ChromEvolOptions::globalParams_);
    std::map<int, bool> mapOfInterSharedModels;
    // removing the previous examined model
    while(paramNumIt != sharedParams.end()){
        int paramNum = paramNumIt->first;
        // how many shared parameters in this particular category
        size_t sharedParamsSize = sharedParams[paramNum].size();
        if (sharedParamsSize == 0){
            auto elemToDel = paramNumIt;
            paramNumIt ++;
            sharedParams.erase(elemToDel);     
            continue;
        }
        // added as the last one
        bool isInterModelShared = false;
        uint firstModel = sharedParams[paramNum][0].first;
        int i = static_cast<int>(sharedParamsSize)-1;
        while (i >=  0){
            uint model = sharedParams[paramNum][static_cast<size_t>(i)].first;
            if (model == numOfShifts){
                sharedParams[paramNum].pop_back();
            }else{
                if (model != firstModel){
                    isInterModelShared = true;
                    mapOfInterSharedModels[paramNum] = isInterModelShared;
                    break;
                }

            }
            i--;
  
        }
        if (sharedParams[paramNum].size() == 0){
            auto elemToDel = paramNumIt;
            paramNumIt ++;
            sharedParams.erase(elemToDel);     
            continue;
                    
        }
        if (paramNum > maxNum){
            maxNum = paramNum;
        }


        paramNumIt ++;
    }
    
    // update the shared parameters according to the appropriate model from which the new one derives.
    auto it = sharedParams.begin();
    auto numOfParamNums = sharedParams.size();
    size_t currNumOfItems = 0;
    while (it != sharedParams.end()){
        auto paramNum = it->first;
        auto sharedModelAndType = sharedParams[paramNum];
        size_t sizeOfSharedParams = sharedModelAndType.size();
        bool isSharedBetweenModels = mapOfInterSharedModels[paramNum];
        bool maxNumUpdated = false;
        for (size_t i = 0; i < sizeOfSharedParams; i++){
            if (sharedModelAndType[i].first == prevShift){
                std::pair<uint, int> modelAndParam;
                modelAndParam.first = numOfShifts;
                modelAndParam.second = sharedModelAndType[i].second;
                if (isSharedBetweenModels){
                // if for example gain2 = gain1, if model 3 derives from model2, gain3 = gain1.
                    sharedParams[paramNum].push_back(modelAndParam);
                    maxNumUpdated = true;

                }else{

                    auto globalParamIt = std::find(globalParamTypes.begin(), globalParamTypes.end(), modelAndParam.second);
                    if (globalParamIt == globalParamTypes.end()){
                        
                        //paramNumToPut = maxNum + 1;
                        if (!maxNumUpdated){
                            maxNum ++;
                            
                        }
                        sharedParams[maxNum].push_back(modelAndParam);
                        maxNumUpdated = true;
                        
                    
                    }else{
                        // a special case for on model
                        if (std::find(globalsAccounted.begin(), globalsAccounted.end(), modelAndParam.second) == globalsAccounted.end()){
                            globalsAccounted.push_back(modelAndParam.second);

                        }
                        maxNumUpdated = true;
                        sharedParams[paramNum].push_back(modelAndParam);
                    }
                

                }
            }

        }       
        it ++;
        currNumOfItems ++;
        if (currNumOfItems == numOfParamNums){
            break;
        }
    }
    if (numOfShifts <= 2){
        for (size_t i = 0; i < globalParamTypes.size(); i++){
            if (std::find(globalsAccounted.begin(), globalsAccounted.end(), globalParamTypes[i]) == globalsAccounted.end()){
                for (size_t j = 1; j <= numOfShifts; j++){
                    std::pair<uint, int> modelAndParam;
                    modelAndParam.first = static_cast<uint>(j);
                    modelAndParam.second = globalParamTypes[i];
                    sharedParams[maxNum + 1].push_back(modelAndParam);
                

                }
                maxNum ++;
            
            }
        }

    }


}
// /*******************************************************************************/
// void LikelihoodUtils::updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>>* paramNameAndType, SingleProcessPhyloLikelihood* tl, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::string suffix){
//     ParameterList substitutionModelParams = tl->getSubstitutionModelParameters();
//     std::vector<std::string> namesAllParams = substitutionModelParams.getParameterNames();
//     updateMapsOfParamTypesAndNames(typeWithParamNames, paramNameAndType, namesAllParams, sharedParams, suffix);


// }

/**********************************************************************************/
// void LikelihoodUtils::updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>>* paramNameAndType, SingleProcessPhyloLikelihood* tl, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams){
//     ParameterList substitutionModelParams = tl->getSubstitutionModelParameters();
//     std::vector<std::string> namesAllParams = substitutionModelParams.getParameterNames();
//     std::map<string, vector<std::pair<uint, int>>> sharedParamsNames;
//     if (sharedParams){
//         LikelihoodUtils::createMapOfSharedParameterNames(*sharedParams, sharedParamsNames);
//     }
//     for (size_t i = 0; i < namesAllParams.size(); i++){
//         uint modelId = LikelihoodUtils::getModelFromParamName(namesAllParams[i]);
//         int type = LikelihoodUtils::getTypeOfParamFromParamName(namesAllParams[i]);
//         //should get the type

        
//         typeWithParamNames[type][modelId].push_back(namesAllParams[i]);
//         if (paramNameAndType){
//             (*paramNameAndType)[namesAllParams[i]] = std::pair<int, uint>(type, modelId);
//         }
//         if (sharedParams){
//             if (sharedParamsNames.find(namesAllParams[i]) != sharedParamsNames.end()){
//                 for (size_t k = 0; k < sharedParamsNames[namesAllParams[i]].size(); k++ ){
//                     uint model = sharedParamsNames[namesAllParams[i]][k].first;
//                     int typeOfSharedParam = sharedParamsNames[namesAllParams[i]][k].second;
//                     typeWithParamNames[typeOfSharedParam][model].push_back(namesAllParams[i]);

//                 }

//             }

//         }

//     }

// }
/**********************************************************************************/
void LikelihoodUtils::createMapOfSharedParameterNames(std::map<int, std::vector<std::pair<uint, int>>> &sharedParams, std::map<string, vector<std::pair<uint, int>>> &sharedParamsNames){
    auto it = sharedParams.begin();
    while (it != sharedParams.end()){
        auto sharedPerParamNum = sharedParams[it->first];
        uint model = sharedPerParamNum[0].first;
        int type = sharedPerParamNum[0].second;
        string paramBasicName = LikelihoodUtils::getStringParamName(type);
        string paramName;
        size_t numOfParams = LikelihoodUtils::getNumberOfParametersPerParamType(type, ChromEvolOptions::rateChangeType_);
        for (size_t k = 0; k < numOfParams; k++){
            if (type == ChromosomeSubstitutionModel::BASENUM){
                paramName = "Chromosome." + paramBasicName +"_"+ std::to_string(model);

            }else{
                paramName = "Chromosome." + paramBasicName +std::to_string(k) + "_"+ std::to_string(model);
            }
            for (size_t i = 1; i < sharedPerParamNum.size(); i++){
                sharedParamsNames[paramName].push_back(sharedPerParamNum[i]);
            }
        }

        it ++;
    }

}

/***********************************************************************************************/
void LikelihoodUtils::getMutableMapOfModelAndNodeIds(std::map<uint, vector<uint>> &mapModelNodesIds, SingleProcessPhyloLikelihood* lik, uint rootId){
    uint numOfModels = static_cast<uint>(lik->getSubstitutionProcess().getNumberOfModels());
    if (rootId){
        mapModelNodesIds[1].push_back(rootId);
    }
    for (uint i = 1; i <= numOfModels; i++){
        auto vectorOfNodes = lik->getSubstitutionProcess().getNodesWithModel(i);
        for (size_t j = 0; j < vectorOfNodes.size(); j++){
            mapModelNodesIds[i].push_back(vectorOfNodes[j]);
        }
    }
}
/****************************************************************/
std::map<uint, pair<int, std::map<int, std::vector<double>>>> LikelihoodUtils::getMapOfParamsForComplexModel(SingleProcessPhyloLikelihood* lik, std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames, uint numOfModels) {
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> heterogeneousModelParams;
    for (uint i = 1; i <= numOfModels; i++){
        heterogeneousModelParams[i] = pair<int, std::map<int, std::vector<double>>>();
    }
    auto it = typeWithParamNames.begin();
    while (it != typeWithParamNames.end()){
        int type = it->first;
        auto modelIt = typeWithParamNames[type].begin();
        while(modelIt != typeWithParamNames[type].end()){
            uint model = modelIt->first;
            if (type == ChromosomeSubstitutionModel::BASENUM){
                heterogeneousModelParams[model].first = static_cast<int>(lik->getParameter(typeWithParamNames[type][model][0]).getValue());
            }else{
                vector<string> paramNames = typeWithParamNames[type][model];
                for (size_t i = 0; i < paramNames.size(); i++){
                    double paramValue = lik->getParameter(paramNames[i]).getValue();
                    heterogeneousModelParams[model].second[type].push_back(paramValue);
                }
                
            }
            modelIt ++;
        }
        it ++;
    }
    auto modelIterator = heterogeneousModelParams.begin();
    // it is important to set base number as ignored if it is not used. Other parameters are manipulated
    // by the function specification.
    while (modelIterator != heterogeneousModelParams.end()){
        uint model = modelIterator->first;
        if ((heterogeneousModelParams[model].second[ChromosomeSubstitutionModel::BASENUMR]).size() == 0){
            heterogeneousModelParams[model].first = IgnoreParam;
        }
        modelIterator ++;
    }
    return heterogeneousModelParams;

}
/***********************************************/
void LikelihoodUtils::setParamsNameInForMultiProcess(std::map<uint, std::map<int, vector<string>>> &mapOfParamsNamesPerModelType, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams){
    auto it = modelParams.begin();
    while(it != modelParams.end()){
        uint model = it->first;
        for (int i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
            vector<string> paramsFullNames;
            auto strParamName = LikelihoodUtils::getStringParamName(i);
            if (i == ChromosomeSubstitutionModel::BASENUM){
                auto fullParamName = "Chromosome." + strParamName + "_"+ std::to_string(model);
                mapOfParamsNamesPerModelType[model][i].push_back(fullParamName);
                continue;
            }
            for (size_t j = 0; j < modelParams[model].second[i].size(); j++){
                auto fullParamName = "Chromosome." + strParamName + std::to_string(j)+ "_"+ std::to_string(model);
                mapOfParamsNamesPerModelType[model][i].push_back(fullParamName);
            }
            
        }
        it ++;
    }

}


/**********************************************************************************************/
void LikelihoodUtils::aliasParametersInSubstitutionProcess(std::map<uint, std::map<int, vector<string>>> &mapOfParamsNamesPerModelType, std::map<int, vector<std::pair<uint, int>>>* updatedSharedParams, std::shared_ptr<NonHomogeneousSubstitutionProcess> process){
    std::map<uint, std::map<int, vector<string>>> mapOfUpdatedParamNames;
    auto paramNumIt = updatedSharedParams->begin();
    while (paramNumIt != updatedSharedParams->end()){
        auto paramNum = paramNumIt->first;
        auto pairsOfModelsAndTypes = (*updatedSharedParams)[paramNum]; // vector<pair<uint, int>>
        if (pairsOfModelsAndTypes.size() < 2){
            paramNumIt ++;
            continue;
        }
        uint firstModel = pairsOfModelsAndTypes[0].first;// vector<pair<uint, int>>[0]-> pair<uint, int>.first->uint
        int type = pairsOfModelsAndTypes[0].second;
        
        auto parametersOfFirstModel = mapOfParamsNamesPerModelType[firstModel][type];
        mapOfUpdatedParamNames[firstModel][type] = parametersOfFirstModel;
        for (size_t i = 1; i < pairsOfModelsAndTypes.size(); i++){ 
            uint modelToAlias = pairsOfModelsAndTypes[i].first;
            int typeToAlias = pairsOfModelsAndTypes[i].second;
            auto paramsToAlias = mapOfParamsNamesPerModelType[modelToAlias][typeToAlias]; 
            for (size_t j = 0; j < paramsToAlias.size(); j++){
                process->aliasParameters(parametersOfFirstModel[j],paramsToAlias[j]);

            }           
        }
        paramNumIt ++;
    }
}


/***********************************************************************************/
SubstitutionProcess* LikelihoodUtils::setChromosomeSubstitutionModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, std::map<int, vector<std::pair<uint,int>>>* updatedSharedParams, bool weightedRootFreqs, vector<std::shared_ptr<ChromosomeSubstitutionModel>>* models, std::shared_ptr<ParametrizablePhyloTree> parTree){
    std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    if (parTree == nullptr){
        parTree = std::make_shared<ParametrizablePhyloTree>(*tree);

    }
    //std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
    std::map<uint, std::map<int, vector<string>>> mapOfParamsNamesPerModelType;
    LikelihoodUtils::setParamsNameInForMultiProcess(mapOfParamsNamesPerModelType, modelParams);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, modelParams[1].second, modelParams[1].first, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
    if (weightedRootFreqs){
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree);

    }else{
        vector <double> rootFreqs = LikelihoodUtils::setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
        std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqs);
        std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies);
    }
    
    // adding models
    for (uint i = 1; i <= numOfModels; i++){
        if (i > 1){
            chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, modelParams[i].second, modelParams[i].first, baseNumberUpperBound[i], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
        }   
        subProSim->addModel(chrModel, mapModelNodesIds[i]);
        if (models){
            models->push_back(chrModel);
        }
    }
    LikelihoodUtils::aliasParametersInSubstitutionProcess(mapOfParamsNamesPerModelType, updatedSharedParams, subProSim);


    SubstitutionProcess* nsubPro= subProSim->clone();
    return nsubPro;

}
/***********************************************************************************/
bool LikelihoodUtils::getIfWeightedRootFreq(){
    if (ChromEvolOptions::fixedFrequenciesFilePath_ == "none"){
        return true;
    }
    return false;


}
/***********************************************************************************/
SubstitutionProcess* LikelihoodUtils::setRandomChromosomeSubstitutionModel(std::shared_ptr<ParametrizablePhyloTree> parTree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, double parsimonyBound, std::map<uint, vector<int>> &fixedParams, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, bool weightedRootFreqs, vector<std::shared_ptr<ChromosomeSubstitutionModel>>* models){
    //std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
    std::map<uint, std::map<int, vector<string>>> mapOfParamsNamesPerModelType;
    LikelihoodUtils::setParamsNameInForMultiProcess(mapOfParamsNamesPerModelType, modelParams);
    std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::shared_ptr<ChromosomeSubstitutionModel>(ChromosomeSubstitutionModel::initRandomModel(alphabet, modelParams[1].first, modelParams[1].second, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_, fixedParams[1], parsimonyBound));
    if (weightedRootFreqs){
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree);

    }else{
        vector <double> rootFreqs = LikelihoodUtils::setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
        std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqs);
        std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies);
    }

    
    // adding models
    for (uint i = 1; i <= numOfModels; i++){
        if (i > 1){
            chrModel = std::shared_ptr<ChromosomeSubstitutionModel>(ChromosomeSubstitutionModel::initRandomModel(alphabet, modelParams[i].first, modelParams[i].second, baseNumberUpperBound[i], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_, fixedParams[i], parsimonyBound));
        }  
        subProSim->addModel(chrModel, mapModelNodesIds[i]);
        if (models){
            models->push_back(chrModel);
        }
    }
    LikelihoodUtils::aliasParametersInSubstitutionProcess(mapOfParamsNamesPerModelType, sharedParams, subProSim);
    SubstitutionProcess* nsubPro= subProSim->clone();
    return nsubPro;
}
/***********************************************************************************/
SubstitutionProcess* LikelihoodUtils::setRandomChromosomeSubstitutionModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, double parsimonyBound, std::map<uint, vector<int>> &fixedParams, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, bool weightedRootFreqs, vector<std::shared_ptr<ChromosomeSubstitutionModel>>* models, std::shared_ptr<ParametrizablePhyloTree> parTree){
    if (parTree == nullptr){
        parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
    }
    SubstitutionProcess* substitutionModel = setRandomChromosomeSubstitutionModel(parTree, vsc, alphabet, baseNumberUpperBound, mapModelNodesIds, modelParams, numOfModels, parsimonyBound, fixedParams, sharedParams, weightedRootFreqs, models);
    return substitutionModel;


}


/***********************************************************************************/
// SubstitutionProcess* LikelihoodUtils::setRandomChromosomeSubstitutionModel(PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, double parsimonyBound, std::map<uint, vector<int>> &fixedParams, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, bool weightedRootFreqs, vector<std::shared_ptr<ChromosomeSubstitutionModel>>* models, std::shared_ptr<ParametrizablePhyloTree> parTree){
//     if (parTree == nullptr){
//         parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
//     }
//     SubstitutionProcess* substitutionModel = setRandomChromosomeSubstitutionModel(parTree, vsc, alphabet, baseNumberUpperBound, mapModelNodesIds, modelParams, numOfModels, parsimonyBound, fixedParams, sharedParams, weightedRootFreqs, models);
//     return substitutionModel;


// }
/***********************************************************************************/
SingleProcessPhyloLikelihood* LikelihoodUtils::setHeterogeneousModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, std::map<int, vector<std::pair<uint, int>>>* updatedSharedParams){
    bool weightedRootFreqs = getIfWeightedRootFreq();
    SubstitutionProcess* nsubPro = setChromosomeSubstitutionModel(tree, vsc, alphabet, baseNumberUpperBound, mapModelNodesIds, modelParams, numOfModels, updatedSharedParams, weightedRootFreqs, 0, nullptr);
    Context* context = new Context();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;

}
/**********************************************************************************************/
SingleProcessPhyloLikelihood* LikelihoodUtils::setRandomHeterogeneousModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet,
    std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, 
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, 
    uint numOfModels, double parsimonyBound, 
    std::map<uint, vector<int>> &fixedParams,
    std::map<int, std::vector<std::pair<uint, int>>>* sharedParams)
{
    bool weightedRootFreqs = getIfWeightedRootFreq();
    SubstitutionProcess* nsubPro = setRandomChromosomeSubstitutionModel(tree, vsc,  alphabet, baseNumberUpperBound,mapModelNodesIds, modelParams, numOfModels, parsimonyBound, fixedParams, sharedParams, weightedRootFreqs, 0, nullptr);
    Context* context = new Context();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;

}
/***************************************************************************************************************************************/
void LikelihoodUtils::separateBetweenModels(JointPhyloLikelihood* lik, std::string &traitModel, std::map<uint, std::vector<std::string>> &paramsPerModel){
  string traitModelPrefix;
  string chromosomeModelPrefix = "Chromosome";
  if (traitModel == "Binary"){
    traitModelPrefix = "TwoParameterBinary";

  }else{
    throw Exception("JointTraitChromosomeLikelihood::separateBetweenModels(): only binary trait model is implemented now!");
  }
  auto parameters = lik->getParameters();
  auto paramsNames = parameters.getParameterNames();
  for (size_t i = 0; i < paramsNames.size(); i++){
    if (paramsNames[i].size() >= traitModelPrefix.size() && paramsNames[i].compare(0, traitModelPrefix.size(), traitModelPrefix) == 0){
      paramsPerModel[1].push_back(paramsNames[i]);

    }else if (paramsNames[i].size() >= chromosomeModelPrefix.size() && paramsNames[i].compare(0, chromosomeModelPrefix.size(), chromosomeModelPrefix) == 0){
      paramsPerModel[2].push_back(paramsNames[i]);
    }
  }
  return;

}
