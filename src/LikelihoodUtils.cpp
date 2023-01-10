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