#include "StochasticMappingUtils.h"


using namespace bpp;
using namespace std;

void StochasticMappingUtils::printResultsForEachMapping(PhyloTree* tree, std::map<uint, std::map<int, double>> &expectationsPerTypeRootToNode, const NonHomogeneousSubstitutionProcess* NonHomoProcess, std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> &rootToLeafTransitions, std::map<uint, std::map<size_t, bool>> &presentMapping, const string &outStMappingRootToLeafPath, size_t mappingIndex){
    ofstream stream;
    std::map<uint, size_t> modelsForBranch = getModelForEachBranch(*tree, *NonHomoProcess);
    stream.open(outStMappingRootToLeafPath);
    
    stream << "NODE";
    for (int type = 0; type < ChromosomeSubstitutionModel::typeOfTransition::NUMTYPES; type++){
        string typeStr = getTypeOfTransitionStr(type);
        stream << "," << typeStr;
    }
    stream << std::endl;
    auto nodes = tree->getAllNodes();
    // expectationsPerTypeRootToLeaf
    for (size_t i = 0; i < nodes.size(); i++){
        uint nodeId = tree->getNodeIndex(nodes[i]);
        if (nodeId == tree->getRootIndex()){
            continue;
        }
        if (tree->isLeaf(nodeId)){
            stream << nodes[i]->getName();
        }else{
            stream << "N" << nodeId;
        }
        
        if (!(presentMapping[nodeId][mappingIndex])){
            writeNanInTable(stream);
        }else{
            if (rootToLeafTransitions.find(nodeId) == rootToLeafTransitions.end()){
                writeZeroInTable(stream);
                if (expectationsPerTypeRootToNode.find(nodeId) == expectationsPerTypeRootToNode.end()){
                    for (int type = 0; type < ChromosomeSubstitutionModel::typeOfTransition::NUMTYPES; type++){
                        expectationsPerTypeRootToNode[nodeId][type] = 0;
                    }

                }

            }else{
                if (rootToLeafTransitions[nodeId].find(mappingIndex) == rootToLeafTransitions[nodeId].end()){
                    writeZeroInTable(stream);
                    if (expectationsPerTypeRootToNode.find(nodeId) == expectationsPerTypeRootToNode.end()){
                        for (int type = 0; type < ChromosomeSubstitutionModel::typeOfTransition::NUMTYPES; type++){
                            expectationsPerTypeRootToNode[nodeId][type] = 0;
                        }
                    }
                }else{
                    auto &transitions = rootToLeafTransitions[nodeId][mappingIndex];
                    auto model = NonHomoProcess->getModel(modelsForBranch[nodeId]);
                    auto chrModel = dynamic_cast<const ChromosomeSubstitutionModel*>(model);
                    std::map<int, double> transitionsPerType = getTypeForEachTransitionPerNode(chrModel, transitions, nodeId);
                    for (int type = 0; type < ChromosomeSubstitutionModel::typeOfTransition::NUMTYPES; type++){
                        stream << "," << transitionsPerType[type];
                        if (expectationsPerTypeRootToNode.find(nodeId) == expectationsPerTypeRootToNode.end()){
                            expectationsPerTypeRootToNode[nodeId][type] = transitionsPerType[type];
                        }else{
                            expectationsPerTypeRootToNode[nodeId][type] += transitionsPerType[type];
                        }
                    }
                    stream << std::endl;

                }
            }
            
        }

    }

    stream.close();

}


/**************************************************************************************/
void StochasticMappingUtils::printStochasticMappingEvolutionaryPath(ChromosomeAlphabet* alphabet, std::shared_ptr<PhyloTree> stmTree, const std::map<uint, std::vector<MutationPath>> &mappings, const std::map<uint, std::vector<size_t>> &ancestralStates, size_t mappingIndex, const string &outPathPerMapping){
    ofstream outFile;
    outFile.open(outPathPerMapping);
    size_t totalNumTransitions = 0;
    vector<shared_ptr<PhyloNode> > nodes = stmTree->getAllNodes();
    size_t nbNodes = nodes.size();
    for (size_t n = 0; n < nbNodes; n++){
        uint nodeId = stmTree->getNodeIndex(nodes[n]);
        if (stmTree->getRootIndex() == nodeId){
            outFile << "N-" + std::to_string(nodeId) << endl;
            size_t rootState = ancestralStates.at(nodeId)[mappingIndex];
            outFile <<"\tThe root state is: "<< ((int)(rootState + alphabet->getMin())) <<endl;
        }else{
            if (stmTree->isLeaf(nodeId)){
                outFile << stmTree->getNode(nodeId)->getName() << endl;
            }else{
                outFile << "N-" + std::to_string(nodeId) <<endl;

            }
            MutationPath mutPath = mappings.at(nodeId)[mappingIndex];
            vector<size_t> states = mutPath.getStates();
            vector<double> times = mutPath.getTimes();
            totalNumTransitions += static_cast<int>(times.size());

            auto edgeIndex =  stmTree->getIncomingEdges(nodeId)[0]; 
            auto fatherIndex = stmTree->getFatherOfEdge(edgeIndex);
            outFile << "Father is: " << "N-" << fatherIndex << std::endl;
            size_t fatherState = ancestralStates.at(fatherIndex)[mappingIndex] + alphabet->getMin();    
            for (size_t i = 0; i < states.size(); i++){
                outFile << "from state: "<< fatherState  <<"\tt = "<< times[i] << " to state = "<< ((int)(states[i]) + alphabet->getMin()) << endl;
                fatherState = ((int)(states[i]) + alphabet->getMin());
            }
            outFile <<"# Number of transitions per branch: "<< times.size() <<endl;   

        }

        outFile <<"*************************************"<<endl;

    }
    outFile <<"Total number of transitions is: "<< totalNumTransitions << endl;
    outFile.close();

}
/**************************************************************************************/
void StochasticMappingUtils::printMappingStretchedBranches(const string &stretchedBranchesPath, std::map<uint, std::map<size_t, std::pair<double, double>>> &originalBrLenVsStretched, const std::shared_ptr<PhyloTree> tree){
    if (originalBrLenVsStretched.empty()){
        return;
    }
    ofstream outFile;
    outFile.open(stretchedBranchesPath);
    outFile << "Node,Mapping,OriginalBranchLength,StretchedBranchLength" << std::endl;
    auto itNode = originalBrLenVsStretched.begin();
    while (itNode != originalBrLenVsStretched.end()){
        auto itMapping = originalBrLenVsStretched[itNode->first].begin();
        while (itMapping != originalBrLenVsStretched[itNode->first].end()){
            if(tree->isLeaf(tree->getNode(itNode->first))){
                outFile << (tree->getNode(itNode->first))->getName() << ",";
            }else{
                outFile << "N" << itNode->first << ",";
            }
            outFile << itMapping->first << ",";
            outFile << (itMapping->second).first << ",";
            outFile << (itMapping->second).second << std::endl;
            itMapping++;
        }
        itNode ++;
    }

    outFile.close();

}
/**************************************************************************************/
void StochasticMappingUtils::fixFailedMappings(PhyloTree* tree, StochasticMapping* stm, std::map<uint, std::map<size_t, std::pair<double, double>>> &originalBrLenVsStretched, size_t numOfFixingMappingIterations){
    auto failedNodesWithMappings = stm->getFailedNodes();
    std::cout << "*** Stochastic mapping: failed nodes before heuristics ***"<< std:: endl;
    auto it = failedNodesWithMappings.begin();
    while (it != failedNodesWithMappings.end()){
        if (tree->isLeaf(it->first)){
            std::cout << (tree->getNode(it->first))->getName() <<": ";
        }else{
            std::cout << "N" << it->first << ": ";
        }
        for (size_t m = 0; m < failedNodesWithMappings[it->first].size(); m++){
            if (m == failedNodesWithMappings[it->first].size()-1){
                std::cout << failedNodesWithMappings[it->first][m] << std::endl;

            }else{
                std::cout << failedNodesWithMappings[it->first][m] <<", ";
            }
            
        }
        it++;
    }
    std::cout << "******End of unrepresented nodes*********" << std::endl;
    vector <uint> failedNodes = getVectorOfMapKeys(failedNodesWithMappings);
    for (size_t i = 0; i < failedNodes.size(); i++){
        auto mappings = failedNodesWithMappings[failedNodes[i]];
        for (size_t j = 0; j < mappings.size(); j++){
            auto branchPtr = tree->getIncomingEdges(tree->getNode(failedNodes[i]))[0];
            double original_branch_length = branchPtr->getLength();
            auto branchLength = original_branch_length;
            for (size_t k = 0; k < MAX_ITER_HEURISTICS; k++){
                auto rateToLeave = stm->getRateToLeaveState(failedNodes[i], mappings[j]);
                if (branchLength * rateToLeave >= 1){
                    branchLength *= BRANCH_MULTIPLIER_FACTOR;
                }else{
                    branchLength *= 1/(rateToLeave*branchLength);
                }

                bool success = stm->tryToReplaceMapping(branchLength, failedNodes[i], mappings[j], numOfFixingMappingIterations);
                if (success){
                    stm->removeFailedNodes(failedNodes[i], mappings[j]);
                    originalBrLenVsStretched[failedNodes[i]][mappings[j]] = std::pair<double, double>(original_branch_length, branchLength);
                    break;
                }
            }           
        }
    }
}
/**************************************************************************************/
void StochasticMappingUtils::printRootToLeaf(PhyloTree* tree, ChromosomeAlphabet* alphabet, std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> &rootToLeafOccurrences, std::map<uint, std::map<size_t, bool>> &presentMapping, size_t numOfMappings, const NonHomogeneousSubstitutionProcess* NonHomoProcess, const string &resultsDir){
  const string outStMappingRootToLeafPath =  resultsDir+"//"+ "stMapping_root_to_leaf.txt";
  const string outStMappingRootToLeafExpPath =  resultsDir+"//"+ "stMapping_root_to_leaf_exp.csv";
  ofstream stream;
  stream.open(outStMappingRootToLeafPath);
  std::map<uint, std::map<std::pair<size_t, size_t>, double>> expectationsFromRootToLeaf;
  stream << "############################" << std::endl;
  stream << "# Root to leaf transitions #" << std::endl;
  stream << "############################" << std::endl;
  std::map<uint, size_t> modelsForBranch = getModelForEachBranch(*tree, *NonHomoProcess);
  auto leaves = tree->getAllLeaves();
  std::map<uint, double> numOfMappingsPerLeaf;
  for (size_t i = 0; i < leaves.size(); i++){
    double numOfAccountedMappings = 0;
    auto leafIndex = tree->getNodeIndex(leaves[i]);
    stream << "*** Leaf name: " << leaves[i]->getName() << std::endl;
    for (size_t j = 0; j < numOfMappings; j ++){
      if (presentMapping[leafIndex][j]){
        stream << "\t# Mapping $" << j << std::endl;
        if (rootToLeafOccurrences[leafIndex].find(j) == rootToLeafOccurrences[leafIndex].end()){
            if (expectationsFromRootToLeaf.find(leafIndex) == expectationsFromRootToLeaf.end()){
                expectationsFromRootToLeaf[leafIndex];
            }
            numOfAccountedMappings ++;
            continue;
        }
        numOfAccountedMappings ++;

        auto &transitions = rootToLeafOccurrences[leafIndex][j];
        auto it = transitions.begin();
        while(it != transitions.end()){
            if (expectationsFromRootToLeaf.find(leafIndex) == expectationsFromRootToLeaf.end()){
                expectationsFromRootToLeaf[leafIndex][it->first] = transitions[it->first];
            }else{
                if (expectationsFromRootToLeaf[leafIndex].find(it->first) == expectationsFromRootToLeaf[leafIndex].end()){
                    expectationsFromRootToLeaf[leafIndex][it->first] = transitions[it->first];
                }else{
                    expectationsFromRootToLeaf[leafIndex][it->first] += transitions[it->first];
                }
            }
            stream <<"\t\t" << alphabet->getMin() + (it->first).first << " -> " << alphabet->getMin()+ (it->first).second << " : " <<  transitions[it->first] << std::endl;
            it ++;
        }

      }
    }
    if (expectationsFromRootToLeaf.find(leafIndex) != expectationsFromRootToLeaf.end()){
       auto &expTransitios = expectationsFromRootToLeaf[leafIndex]; 
       auto transitionsIt = expTransitios.begin();
       while(transitionsIt != expTransitios.end()){
           expectationsFromRootToLeaf[leafIndex][transitionsIt->first] /= numOfAccountedMappings;
           transitionsIt ++;
        }
    }
    numOfMappingsPerLeaf[leafIndex] = numOfAccountedMappings;
  }
  stream << "###############################################" << std::endl;
  stream << "# Root to leaf  transitions expectations      #" << std::endl;
  stream << "###############################################" << std::endl;
  for (size_t i = 0; i < leaves.size(); i++){
    auto leafIndex = tree->getNodeIndex(leaves[i]);
    if (expectationsFromRootToLeaf.find(leafIndex) == expectationsFromRootToLeaf.end()){
        continue;
    }
    stream << "*** Leaf name: " << leaves[i]->getName() << std::endl;
    auto &transitionsPerLeaf = expectationsFromRootToLeaf[leafIndex];
    auto itTransitionsExpRootLeaf = transitionsPerLeaf.begin();
    while (itTransitionsExpRootLeaf != transitionsPerLeaf.end()){
        stream << "\t" << alphabet->getMin() + (itTransitionsExpRootLeaf->first).first << " -> " << alphabet->getMin()+ (itTransitionsExpRootLeaf->first).second << " : " <<  transitionsPerLeaf[itTransitionsExpRootLeaf->first] << std::endl;
        itTransitionsExpRootLeaf ++;
    }
  }
  stream.close();
  std::map<uint, std::map<int, double>> expectationsPerTypeRootToLeaf;
  
  for (size_t i = 0; i < numOfMappings; i++){
      const string outStMappingRootToLeafCSVPath =  resultsDir+"//"+ "stMapping_mapping_" + std::to_string(i) + ".csv";
      printResultsForEachMapping(tree, expectationsPerTypeRootToLeaf, NonHomoProcess, rootToLeafOccurrences, presentMapping, outStMappingRootToLeafCSVPath, i);
  }
  ofstream stream_exp;
  stream_exp.open(outStMappingRootToLeafExpPath);
  stream_exp << "NODE";
  for (int type = 0; type < ChromosomeSubstitutionModel::typeOfTransition::NUMTYPES; type++){
    string typeStr = getTypeOfTransitionStr(type);
    stream_exp << "," << typeStr;
  }
  stream_exp << std::endl;
  for(size_t i = 0; i < leaves.size(); i++){
    stream_exp << leaves[i]->getName();
    auto leafIndex = tree->getNodeIndex(leaves[i]);
    if (expectationsPerTypeRootToLeaf.find(leafIndex) == expectationsPerTypeRootToLeaf.end()){
        writeNanInTable(stream_exp);
    }else{

        for (int type = 0; type < ChromosomeSubstitutionModel::typeOfTransition::NUMTYPES; type++){
          expectationsPerTypeRootToLeaf[leafIndex][type] /= numOfMappingsPerLeaf[leafIndex];
          stream_exp << "," << expectationsPerTypeRootToLeaf[leafIndex][type];

        }
        stream_exp << std::endl;

    }

  }

  stream_exp.close();
}

/*************************************************************************************/
void StochasticMappingUtils::writeNanInTable(ofstream &stream){
    for (int type = 0; type < ChromosomeSubstitutionModel::typeOfTransition::NUMTYPES; type++){
        stream << ",";
    }
    stream << std::endl;

}

/*************************************************************************************/
void StochasticMappingUtils::writeZeroInTable(ofstream &stream){
    for (int type = 0; type < ChromosomeSubstitutionModel::typeOfTransition::NUMTYPES; type++){
        stream << ",0";
    }
    stream << std::endl;
}
/**************************************************************************************/
std::map<uint, size_t> StochasticMappingUtils::getModelForEachBranch(PhyloTree &tree, const NonHomogeneousSubstitutionProcess &NonHomoModel){
    std::map<uint, size_t> modelPerBranch;
    auto rootId = tree.getRootIndex();
    auto sons = tree.getSons(rootId);
    for (size_t i = 0; i < sons.size(); i++){
        getModeForSons(tree, rootId, sons[i], &NonHomoModel, modelPerBranch);

    }
    return modelPerBranch;
}

/*************************************************************************************/
void StochasticMappingUtils::getModeForSons(PhyloTree &tree, uint fatherId, uint nodeId, const NonHomogeneousSubstitutionProcess* NonHomoModel, std::map<uint, size_t> &modelPerNode){
    if (fatherId == tree.getRootIndex()){
        modelPerNode[nodeId] = 1;
    }else{
        modelPerNode[nodeId] = NonHomoModel->getModelNumberForNode(fatherId);
    }
    if (tree.isLeaf(nodeId)){
        return;
    }
    auto sons = tree.getSons(nodeId);
    for (size_t i = 0; i < sons.size(); i++){
        getModeForSons(tree, nodeId, sons[i], NonHomoModel, modelPerNode);
    }

}
/*************************************************************************************/

std::string StochasticMappingUtils::getTypeOfTransitionStr(int transitionType){
    string nameOfTransition;
    if (transitionType == ChromosomeSubstitutionModel::GAIN_T){
        nameOfTransition = "GAIN";
    }else if (transitionType == ChromosomeSubstitutionModel::LOSS_T){
        nameOfTransition = "LOSS";
    }else if (transitionType == ChromosomeSubstitutionModel::DUPL_T){
        nameOfTransition =  "DUPLICATION";
    }else if (transitionType == ChromosomeSubstitutionModel::DEMIDUPL_T){
        nameOfTransition = "DEMI-DUPLICATION";
    }else if (transitionType == ChromosomeSubstitutionModel::BASENUM_T){
        nameOfTransition = "BASE-NUMBER";
    }else if (transitionType == ChromosomeSubstitutionModel::MAXCHR_T){
        nameOfTransition = "TOMAX";
    }else{
        throw Exception("getTypeOfTransitionStr(): No such transition!!");
    }
    return nameOfTransition;
    
}
/**************************************************************/
std::map<int, double> StochasticMappingUtils::getTypeForEachTransitionPerNode(const ChromosomeSubstitutionModel* chrModel, std::map<pair<size_t, size_t>, double> &transitionsPerNode, uint nodeId){
    std::map<int, double> expectationsPerType;
    auto itTransitions = transitionsPerNode.begin();
    while(itTransitions != transitionsPerNode.end()){
        std::vector<double> probabilities;
        int startState = static_cast<int>((itTransitions->first).first);
        int endState = static_cast<int>((itTransitions->first).second);
        double expectation = transitionsPerNode[itTransitions->first];
        bool legalMove = getProbabilitiesPerType(probabilities, startState, endState, chrModel);
        if (!legalMove){
            itTransitions ++;
            continue;
        }
        for (int i = 0; i < ChromosomeSubstitutionModel::typeOfTransition::NUMTYPES; i++){
            if (expectationsPerType.find(i) == expectationsPerType.end()){
                expectationsPerType[i] = 0;
            }
            expectationsPerType[i] += (probabilities[i] * expectation);
        }
        itTransitions ++;
    }
    return expectationsPerType;

}
/**************************************************************/
bool StochasticMappingUtils::getProbabilitiesPerType(vector<double> &probabilities, int startStateIndex, int endStateIndex, const ChromosomeSubstitutionModel* model){
    // convert from state index to real chromsome number
    bool legalMove = false;
    probabilities.resize(ChromosomeSubstitutionModel::NUMTYPES);
    std::fill(probabilities.begin(), probabilities.end(), 0);
    const ChromosomeAlphabet* alphabet = dynamic_cast<const ChromosomeAlphabet*>(model->getAlphabet());
    int chrStart = startStateIndex + alphabet->getMin();
    int chrEnd = endStateIndex + alphabet->getMin();
    double sumOfRates = 0; //for normalization of weights
    //gain
    if (chrStart + 1 == chrEnd){
        if (!(model->isIgnoredGain())){
            if (chrStart > 3){
                probabilities[(size_t)(ChromosomeSubstitutionModel::GAIN_T)] = 1;
                return true;
            }else{
                double gainRate = model->getGain()->getRate(chrStart);
                probabilities[(size_t)(ChromosomeSubstitutionModel::GAIN_T)] = gainRate;
                sumOfRates += gainRate;
                legalMove = true;
            }

        }
    }
    // loss
    if (chrStart - 1 == chrEnd){
        if (!(model->isIgnoredLoss())){
            probabilities[(size_t)(ChromosomeSubstitutionModel::LOSS_T)] = 1;
            return true;
        }
    }
    //baseNumber transitions
    int baseNumber = model->getBaseNumber();
    if (baseNumber != IgnoreParam){
        if (chrEnd > chrStart){
            if (((chrEnd - chrStart) % baseNumber == 0) && ((chrEnd - chrStart) <= (int)(model->getMaxChrRange()))){
                legalMove = true;
                probabilities[(size_t)(ChromosomeSubstitutionModel::BASENUM_T)] = model->getBaseNumR()->getRate(chrStart);
                sumOfRates += model->getBaseNumR()->getRate(chrStart);                              
            }
        }        
    }
    //duplication
    if (chrEnd == 2 * chrStart){
        if (!(model->isIgnoredDupl())){
            legalMove = true;
            double duplRate = model->getDupl()->getRate(chrStart);
            probabilities[(size_t)(ChromosomeSubstitutionModel::DUPL_T)] = duplRate;
            sumOfRates += duplRate;
        }

    }
    //demi-duplication
    if (!(model->isIgnoredDemiDupl())){
        if (chrStart % 2 == 0){
            if (chrEnd == chrStart * 1.5){
                legalMove = true;
                probabilities[(size_t)(ChromosomeSubstitutionModel::DEMIDUPL_T)] = model->getDemiDupl()->getRate(chrStart);
                sumOfRates += model->getDemiDupl()->getRate(chrStart);
            }
        }else{
            if ((chrEnd == (int)ceil(chrStart * 1.5)) || (chrEnd == (int)floor(chrStart * 1.5))){
                legalMove = true;
                double demiDupRate;
                if (chrStart == 1){
                    demiDupRate =  model->getDemiDupl()->getRate(chrStart);
                }else{
                    demiDupRate = model->getDemiDupl()->getRate(chrStart)/2;
                }
                probabilities[(size_t)(ChromosomeSubstitutionModel::DEMIDUPL_T)] = demiDupRate;
                sumOfRates += demiDupRate;
            }
        }

    }
    if (chrEnd  == (int)(alphabet->getMax())){

        
        double toMaxRate = model->Qij(startStateIndex, endStateIndex)-sumOfRates;
        probabilities[(size_t)(ChromosomeSubstitutionModel::MAXCHR_T)] = toMaxRate;
        sumOfRates += toMaxRate;
        if (sumOfRates == 0){
            legalMove = false;
        }else{
            legalMove = true;

        }
    }
    //if nothing fits
    if(!legalMove){
        //throw Exception ("ERROR: ComputeChromosomeTransitionsExp::getProbabilitiesPerType(): Illegal transition!");
        return legalMove;
    }
    for (size_t i = 0; i < probabilities.size(); i++){
        probabilities[i] /= sumOfRates;
    }
    return true;

}

/**************************************************************************************/
vector <uint> StochasticMappingUtils::getVectorOfMapKeys(std::map<uint, vector<size_t>> &mapOfVectors){
    auto it = mapOfVectors.begin();
    vector <uint> vectorOfKeys;
    while (it != mapOfVectors.end()){
        vectorOfKeys.push_back(it->first);
        it ++;
    }
    return vectorOfKeys;
}
