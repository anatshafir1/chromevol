#include "ComputeChromosomeTransitionsExp.h"
using namespace bpp;

/****************************************************************/
void ComputeChromosomeTransitionsExp::computeExpectationOfChangePerBranch(uint nodeId, VVdouble &jointProbFatherNode, int jumpType){

    double expectation = 0;
    auto it = branchTransitionsExp_[nodeId].begin();
    while (it != branchTransitionsExp_[nodeId].end()){
        int fatherState = it->first.first;
        int sonState = it->first.second;
        expectation += jointProbFatherNode[(size_t)sonState][(size_t)fatherState] * getExpectation(nodeId, fatherState, sonState, jumpType);
        it ++;
    }

    expNumOfChangesPerBranch_[nodeId][jumpType] += expectation;
    expNumOfChanges_[jumpType] += expectation;
    return;

}


// /**********************************************************************************/
void ComputeChromosomeTransitionsExp::computeExpPerTypeHeuristics(map <uint, vector<pair<int, int>>>& nonAccountedForBranchesFromFirstRun){
    map <uint, vector<pair<int, int>>>::iterator it = nonAccountedForBranchesFromFirstRun.begin();
    while (it != nonAccountedForBranchesFromFirstRun.end()){
        uint nodeId = it->first;
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
            auto itExistingTerminals =  branchTransitionsExp_[nodeId].begin();
            while(itExistingTerminals != branchTransitionsExp_[nodeId].end()){
                auto accountedTerminal = std::find(nonAccountedForBranchesFromFirstRun[nodeId].begin(), nonAccountedForBranchesFromFirstRun[nodeId].end(), itExistingTerminals->first);
                if (accountedTerminal != nonAccountedForBranchesFromFirstRun[nodeId].end()){
                    // it is an accounted terminal
                    itExistingTerminals++;
                    continue;
                }else{
                    size_t father = itExistingTerminals->first.first;
                    size_t son = itExistingTerminals->first.second;          
                    double expectation = jointProbabilitiesFatherSon_[nodeId][0][son][father] * getExpectation(nodeId, (int)father, (int)son, i);
                    expNumOfChangesPerBranch_[nodeId][i] += expectation;
                    expNumOfChanges_[i] += expectation;
                    itExistingTerminals ++;
                }
            }           
        }
        it ++;
    }



}

// /**********************************************************************************/
void ComputeChromosomeTransitionsExp::computeExpectationPerType(){
    //init();
    vector<shared_ptr<PhyloNode> > nodes = tree_->getAllNodes();
    size_t nbNodes = nodes.size();
    for (size_t n = 0; n < nbNodes; n++){
        uint nodeId = tree_->getNodeIndex(nodes[n]);
        if (tree_->getRootIndex() == nodeId){
            continue;
        }
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
            computeExpectationOfChangePerBranch(nodeId, jointProbabilitiesFatherSon_[nodeId][0], i);

        }
    }

}
/**************************************************************************************/
void ComputeChromosomeTransitionsExp::setNodesOnPathMap(std::map<uint,std::vector<uint>>* nodesOnPath){
    nodesOnPath = new std::map<uint,std::vector<uint>>();
    uint rootId = tree_->getRootIndex();
    auto nodes = tree_->getAllNodes();
    for (size_t i = 0; i < nodes.size(); i++){
        auto name = nodes[i]->getName();
        if (rootId == tree_->getNodeIndex(nodes[i])){
            (*nodesOnPath)[rootId].push_back(rootId);
            continue;
        }
        if (name.find("dummy") != std::string::npos){
            continue;
        }
        uint nodeId = tree_->getNodeIndex(nodes[i]);
        (*nodesOnPath)[nodeId].push_back(nodeId);
        auto fatherNode = tree_->getFatherOfNode(nodes[i]);
        while((fatherNode->getName()).find("dummy") != std::string::npos){
            (*nodesOnPath)[nodeId].push_back(tree_->getNodeIndex(fatherNode));
            fatherNode = tree_->getFatherOfNode(fatherNode);
        }
    }
}

// /*********************************************************************************/
double ComputeChromosomeTransitionsExp::getSumOfAllDummyBranchesOnPath(std::map<uint,std::vector<uint>>* nodesOnPathMap, uint nodeId, int type){
    auto &nodesOnPath = (*nodesOnPathMap)[nodeId];
    double sumOfBranchesOnPath = 0;
    for (size_t i = 0; i < nodesOnPath.size(); i++){
        sumOfBranchesOnPath += expNumOfChangesPerBranch_[nodeId][type];
    }
    return sumOfBranchesOnPath;

}
/*********************************************************************************/
void ComputeChromosomeTransitionsExp::printResults(const string path, bool jointTraitModel) {
    ofstream outFile;
    if (path == "none"){
        throw Exception("ERROR!!! ComputeChromosomeTransitionsExp::printResults(): not provided file path!\n");
    }
    std::map<uint,std::vector<uint>>* nodesOnPath = 0;
    if (jointTraitModel){
        setNodesOnPathMap(nodesOnPath);

    }
    outFile.open(path);
    vector<string> typeNames;
    vector<shared_ptr<PhyloNode> > nodes = tree_->getAllNodes();
    size_t nbNodes = nodes.size();
    for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
        vector<string> statementes;
        string typeOfRate;
        if (i == ChromosomeSubstitutionModel::GAIN_T){
            typeOfRate = "GAIN";
            statementes.push_back("#Nodes with GAIN events with expectation above ");
            statementes.push_back("Total number of gain events: ");
            typeNames.push_back(typeOfRate);
        }else if (i == ChromosomeSubstitutionModel::LOSS_T){
            typeOfRate = "LOSS";
            statementes.push_back("#Nodes with LOSS events with expectation above ");
            statementes.push_back("Total number of loss events: ");
            typeNames.push_back(typeOfRate);
        }else if (i == ChromosomeSubstitutionModel::DUPL_T){
            typeOfRate = "DUPLICATION";
            statementes.push_back("#Nodes with duplication events with expectation above ");
            statementes.push_back("Total number of duplication events: ");
            typeNames.push_back(typeOfRate);
        }else if (i == ChromosomeSubstitutionModel::DEMIDUPL_T){
            typeOfRate = "DEMI-DUPLICATION";
            statementes.push_back("#Nodes with demi-duplication events with expectation above ");
            statementes.push_back("Total number of demi-duplications events: ");
            typeNames.push_back(typeOfRate);
        }else if (i == ChromosomeSubstitutionModel::BASENUM_T){
            typeOfRate = "BASE-NUMBER";
            statementes.push_back("#Nodes with transitions in base number events with expectation above ");
            statementes.push_back("Total number of transitions in base number: ");
            typeNames.push_back(typeOfRate);
        }else{
            typeOfRate = "TOMAX";
            typeNames.push_back(typeOfRate);
            continue;
        }
        outFile << statementes[0] <<THRESHOLD_EXP<<endl;
        
        for (size_t n = 0; n < nbNodes; n++){
            if (jointTraitModel){
                if ((nodes[n]->getName()).find("dummy") !=  std::string::npos){
                    continue;
                }
            }
            uint nodeId = tree_->getNodeIndex(nodes[n]);
            double expectedNumOfChangesPerBranch;
            if (jointTraitModel){
                expectedNumOfChangesPerBranch = getSumOfAllDummyBranchesOnPath(nodesOnPath, nodeId, i);
            }else{
                expectedNumOfChangesPerBranch = expNumOfChangesPerBranch_[nodeId][i];
            }
            if (expectedNumOfChangesPerBranch > THRESHOLD_EXP){
                string nodeName;
                if (tree_->getRootIndex() == nodeId){
                    continue;
                }
                if (tree_->isLeaf(nodeId)){
                    nodeName = (tree_->getNode(nodeId))->getName();

                }else{
                    if (jointTraitModel){
                        nodeName = (tree_->getNode(nodeId))->getName();

                    }else{
                        nodeName = "N" + std::to_string(nodeId);

                    }
                    
                }
                outFile << nodeName <<": "<<expectedNumOfChangesPerBranch<<endl; 
            }
        }
        outFile << statementes[1] <<  expNumOfChanges_[i] <<endl;
        outFile <<"#+++++++++++++++++++++++++++++\n\n";
    }
    outFile << "EVENTS NOT ACCOUNTED FOR IN THE SIMULATIONS: "<< endl;
    for (size_t n = 0; n < nbNodes; n++){
        uint nodeId = tree_->getNodeIndex(nodes[n]);
        if (nodeId == tree_->getRootIndex()){
            continue;
        }
        if (jointTraitModel){
            if ((nodes[n]->getName()).find("dummy") !=  std::string::npos){
                continue;
            }
        }
        double cumulativeProb;
        if (jointTraitModel){
            auto &vectorOfNodesOnPath = (*nodesOnPath)[nodeId];
            string currentNodeName = nodes[n]->getName();
            for (size_t m = 0; m < vectorOfNodesOnPath.size(); m++){
                cumulativeProb = getCumulativeProbability(vectorOfNodesOnPath[m], 0);
                if (cumulativeProb < THRESHOLD_HEURISTIC){
                    auto nodeNameOnPath = tree_->getNode(vectorOfNodesOnPath[m])->getName();
                    outFile << currentNodeName <<"-> "<< nodeNameOnPath <<" : "<< cumulativeProb <<endl;
                }
            }


        }else{
            cumulativeProb = getCumulativeProbability(nodeId, 0);
            if (cumulativeProb < THRESHOLD_HEURISTIC){
                string nodeName;
                if (tree_->isLeaf(nodeId)){
                    nodeName = (tree_->getNode(nodeId))->getName();

                }else{
                    nodeName = "N" + std::to_string(nodeId);
                }
                outFile << nodeName <<": "<< cumulativeProb <<endl;

            }
        }


    }
    outFile <<"#+++++++++++++++++++++++++++++\n\n";
    outFile << "#ALL EVENTS EXPECTATIONS PER NODE"<<endl;
    outFile <<"NODE\t";
    for (size_t i = 0; i < typeNames.size(); i++){
        if (i == typeNames.size()-1){
           outFile <<typeNames[i] << endl;
           continue; 
        }
        outFile <<typeNames[i] <<"\t";
    }
    for (size_t n = 0; n < nbNodes; n++){
        string nodeName;
        uint nodeId = tree_->getNodeIndex(nodes[n]);
        if (jointTraitModel){
            if ((nodes[n]->getName()).find("dummy") != std::string::npos){
                continue;
            }
        }
        if (tree_->getRootIndex() == nodeId){
            continue;
        }
        if (tree_->isLeaf(nodeId)){
            nodeName = (tree_->getNode(nodeId))->getName();
        }else{
            if (jointTraitModel){
                nodeName = (tree_->getNode(nodeId))->getName();
            }else{
                nodeName = "N" + std::to_string(nodeId);

            }
            
        }
        outFile << nodeName <<"\t";
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
            double expectedNumOfChangesPerBranch;
            if (jointTraitModel){
                expectedNumOfChangesPerBranch = getSumOfAllDummyBranchesOnPath(nodesOnPath, nodeId, i);

            }else{
                expectedNumOfChangesPerBranch = expNumOfChangesPerBranch_[nodeId][i];
            }
            
            (i == ChromosomeSubstitutionModel::NUMTYPES-1) ? (outFile << expectedNumOfChangesPerBranch <<endl) : (outFile <<expectedNumOfChangesPerBranch <<"\t");
        }
    }
    outFile << "#+++++++++++++++++++++++++++++\n\n"<<endl;
    //get expected number of changes from root to tip
    outFile << "#EXPECTED NUMBER OF EVENTS FROM ROOT TO LEAF"<<endl;
    outFile <<"NODE\t";
    for (size_t i = 0; i < typeNames.size(); i++){
        if (i == typeNames.size()-1){
           outFile <<typeNames[i] << endl;
           continue; 
        }
        outFile <<typeNames[i] <<"\t";
    }
    auto leaves = tree_->getAllLeaves();
    for (size_t j = 0; j < leaves.size(); j++){
        string leafName = leaves[j]->getName();
        outFile << leafName <<"\t";
        for (int k = 0; k < ChromosomeSubstitutionModel::NUMTYPES; k++){
            double expectedRootToTip = 0;
            uint currNodeId = tree_->getNodeIndex(leaves[j]);
            while (tree_->getRootIndex() != currNodeId){
                expectedRootToTip += expNumOfChangesPerBranch_[currNodeId][k];
                auto currNodeRaw = tree_->getFatherOfNode(tree_->getNode(currNodeId));
                currNodeId = tree_->getNodeIndex(currNodeRaw);
            }
            if (k == ChromosomeSubstitutionModel::NUMTYPES-1){
                outFile << expectedRootToTip <<endl;

            }else{
                outFile << expectedRootToTip <<"\t";
            }
                
        }

    }

    outFile <<"#+++++++++++++++++++++++++++++\n\n";
    outFile <<"#TOTAL EXPECTATIONS:\n";
    // print total expectations
    for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
        outFile << typeNames[i] <<": "<< expNumOfChanges_[i] <<endl;

    }
    if (nodesOnPath){
        delete nodesOnPath;
    }
    outFile.close();
    

}
// /*****************************************************************************************/
PhyloTree* ComputeChromosomeTransitionsExp::getResultTree(const PhyloTree* originalTree){
    std::map<uint,std::vector<uint>>* nodesOnPath = 0;
    if (originalTree){
        setNodesOnPathMap(nodesOnPath);

    }

   PhyloTree* printTree;
   if (originalTree){
    printTree = originalTree->clone();
   }else{
    printTree = tree_->clone();

   }
   vector<shared_ptr<PhyloNode> > nodes = printTree->getAllNodes();
   size_t nbNodes = nodes.size();
    //string branchProp = "expectation";
    
    for (size_t n = 0; n < nbNodes; n++){
        string nodeName;
        uint nodeId = printTree->getNodeIndex(nodes[n]);
        if (printTree->getRootIndex() == nodeId){
            nodeName = "N" + std::to_string(nodeId);
            printTree->getNode(nodeId)->setName(nodeName);
            continue;
        }else if (printTree->isLeaf(nodeId)){
            nodeName = (printTree->getNode(nodeId))->getName();
            
        }else{
            nodeName = "N" + std::to_string(nodeId);

        }
        string expected = "[";
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
            double expectedNumOfChangesPerBranch;
            if (i == ChromosomeSubstitutionModel::NUMTYPES - 1){
                if (originalTree){
                    expectedNumOfChangesPerBranch = getSumOfAllDummyBranchesOnPath(nodesOnPath, nodeId, i);
                }else{
                    expectedNumOfChangesPerBranch = expNumOfChangesPerBranch_[nodeId][i];

                }
                
                expected = expected + to_string(expectedNumOfChangesPerBranch) + "]";
            }else{
                expected = expected + to_string(expectedNumOfChangesPerBranch)+ "\\";
            }            

        }
        nodeName = nodeName + expected;
        printTree->getNode(nodeId)->setName(nodeName);

        
        
    }
    if (nodesOnPath){
        delete nodesOnPath;
    }

    return printTree;

}
// //****************************************************************************************/
void ComputeChromosomeTransitionsExp::init(){
    size_t numOfModels = model_->getNumberOfModels();
    uint rootId = tree_->getRootIndex();
    vector<uint> rootSons = tree_->getSons(rootId);
    for (size_t m = 0; m < numOfModels; m++){
        std::vector<Branch> branchesPerModel;
        vector<std::shared_ptr<PhyloNode>> modelNodes; 
        std::vector<uint> nodeIds = model_->getNodesWithModel(m+1);
        PhyloTree* tree = tree_->clone();
        for (size_t j = 0; j < nodeIds.size(); j++){
            modelNodes.push_back(tree->getNode(nodeIds[j]));
            
        }
        //auto mrca = tree->MRCA(modelNodes);
        auto mrca = TreeUtils::getMRCA(tree, modelNodes);
        if (tree->getNodeIndex(mrca) == tree->getRootIndex()){
            for (size_t k = 0; k < rootSons.size(); k++){
                auto node = tree_->getNode(rootSons[k]);
                auto branchPtr = tree_->getIncomingEdges(node)[0];
                PhyloBranch branch = *branchPtr;
                uint branchIndex = tree_->getEdgeIndex(branchPtr);
                Branch edgeInfo(branchIndex, branch);
                branchesPerModel.push_back(edgeInfo);
                branchTransitionsExp_[rootSons[k]] = std::map<pair<int, int>, pair<int, Vdouble>>();

                for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
                    expNumOfChangesPerBranch_[rootSons[k]][i] = 0;
                }

            }
        }
        delete tree;
        
        for (size_t n = 0; n < nodeIds.size(); n++){
            uint nodeId = nodeIds[n];
            if (nodeId == tree_->getRootIndex()){
                continue;
            }
            if (tree_->isLeaf(nodeId)){
                continue;
            }
            auto node = tree_->getNode(nodeId);
            auto branchPtrs = tree_->getOutgoingEdges(node);
            for (size_t k = 0; k < branchPtrs.size(); k++){
                PhyloBranch branch = *(branchPtrs[k]);
                uint branchIndex = tree_->getEdgeIndex(branchPtrs[k]);
                uint sonId = tree_->getSon(branchIndex);
                Branch edgeInfo(branchIndex, branch);
                branchesPerModel.push_back(edgeInfo);
                
                //std::map<pair<int, int>, pair<int, Vdouble>> terminalsMap();
                branchTransitionsExp_[sonId] = std::map<pair<int, int>, pair<int, Vdouble>>();

                for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
                    expNumOfChangesPerBranch_[sonId][i] = 0;
                }

            }


        }
        branchOrder_.push_back(branchesPerModel);
        sort(branchOrder_[m].begin(), branchOrder_[m].end(), compareBranches);
    }


    for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
        expNumOfChanges_[i] = 0;
    }
    
}

/************************************************************************************/
bool ComputeChromosomeTransitionsExp::compareBranches(Branch& edge1, Branch& edge2){
    return (edge1.second.getLength() < edge2.second.getLength());
}
/************************************************************************************/

void ComputeChromosomeTransitionsExp::runSimulations(int numOfSimulations){
    std::cout << "Initializing expectation instance ..." << std::endl;
    init();
    std::cout << "Run iterations (simulations) ..." << std::endl;
    for (size_t k = 0; k < branchOrder_.size(); k++){
        for (size_t i = 0; i < alphabet_->getSize(); i++){
            for (size_t j = 0; j < (size_t)numOfSimulations; j++){
                runIteration((int)i, k);
            }
        }

    }
    std::cout << "Compute expectations based on simulations ..." << std::endl;
    computeExpectationAndPosterior();
}
// /*************************************************************************************/
void ComputeChromosomeTransitionsExp::runIteration(int beginState, size_t modelIndex, map <uint, vector<pair<int,int>>>* unAccountedNodesAndTerminals){
    auto model = std::dynamic_pointer_cast<const ChromosomeSubstitutionModel>(model_->getModel(modelIndex+1));
    double totalTimeTillJump = 0;
    double maxBranch = branchOrder_[modelIndex][branchOrder_[modelIndex].size()-1].second.getLength();
    int currentState = beginState;
    vector <std::pair<int,int>> jumpsUntilNow;
    // getFatherOfEdge
    int indexOfSmallestNotUpdatedBranch = 0;
    while (totalTimeTillJump < maxBranch){
        double averageWaitingTime = 1 / (-1 * (model->Qij(currentState, currentState)));//waitingTimes_[currentState];
        totalTimeTillJump +=  RandomTools::randExponential(averageWaitingTime);
        for (size_t i = indexOfSmallestNotUpdatedBranch; i < branchOrder_[modelIndex].size(); i++){
            double branchLength = branchOrder_[modelIndex][i].second.getLength();
            if (branchLength > totalTimeTillJump){
                indexOfSmallestNotUpdatedBranch = (int)i;
                break;
            }
            std::pair<int, int> ancestralTerminals;
            ancestralTerminals.first = beginState;
            ancestralTerminals.second = currentState;
            uint branchNodeId = tree_->getSon(branchOrder_[modelIndex][i].first);
            if (unAccountedNodesAndTerminals){
                // unaccountedNodeAndTerminal contain the unaccounted node with their already accounted terminal
                // Therefore, if these terminals are found, we continue       
                if (find((*unAccountedNodesAndTerminals)[branchNodeId].begin(), (*unAccountedNodesAndTerminals)[branchNodeId].end(), ancestralTerminals) != (*unAccountedNodesAndTerminals)[branchNodeId].end()){
                    continue;
                }

            }

            auto terminalsSearch = branchTransitionsExp_[branchNodeId].find(ancestralTerminals);
            if (terminalsSearch == branchTransitionsExp_[branchNodeId].end()){
                branchTransitionsExp_[branchNodeId][ancestralTerminals].first = 0;
                for (size_t k = 0; k < ChromosomeSubstitutionModel::NUMTYPES; k ++){
                    (branchTransitionsExp_[branchNodeId][ancestralTerminals].second).push_back(0);
                }

            }
            branchTransitionsExp_[branchNodeId][ancestralTerminals].first += 1;
            for (size_t j = 0; j < jumpsUntilNow.size(); j++){
                updateExpectationsPerBranch(branchNodeId, ancestralTerminals, jumpsUntilNow[j]);

            }
        }

        int nextState = getRandomState(currentState, model);
        if ((nextState == (int)alphabet_->getSize()-1) && (!isMaxStateValid(currentState, model))){
            break;
        }
        std::pair<int,int> combOcurredStates;
        combOcurredStates.first = currentState;
        combOcurredStates.second = nextState;
        jumpsUntilNow.push_back(combOcurredStates);
        updateMapOfJumps(currentState, nextState, model);
        currentState = nextState;
    }

}
// /********************************************************************************/
bool ComputeChromosomeTransitionsExp::isMaxStateValid(int prevState, std::shared_ptr<const ChromosomeSubstitutionModel> model) const{
    int valid = false;
    int maxState = model->getMax();
    int initState = model->getMin() + prevState;
    // gain
    if (maxState == initState + 1){
        valid = true;
        return valid;
    }
    // dupl
    if (!(model->isIgnoredDupl())){
        if (maxState == 2 * initState){
            valid = true;
            return valid;
        }

    }
    // demi dupl
    if (!(model->isIgnoredDemiDupl())){
        if (initState % 2 == 0){
            if ((int)(initState * 1.5) == maxState){
                valid = true;
            }
        }else{
            if ((maxState == (int)ceil(initState * 1.5)) || (maxState == (int)floor(initState * 1.5))){
                valid = true;
            }
        }
        if (valid){
            return valid;
        }
    }
    // base number
    int baseNumber = model->getBaseNumber();
    if (baseNumber != IgnoreParam){
        if (maxState > initState){
            if (((maxState - initState) % baseNumber == 0) && ((maxState - initState) <= (int)(model->getMaxChrRange()))){
                valid = true;
                return valid;
            }
        }        
    }
    return valid;  
    
}
// /********************************************************************************/
void ComputeChromosomeTransitionsExp::updateMapOfJumps(int startState, int endState, std::shared_ptr<const ChromosomeSubstitutionModel> model){
    bool legalMove = false;
    pair <int, int> jumpStates;
    jumpStates.first = startState;
    jumpStates.second = endState;

    // convert from state index to real chromsome number
    int chrStart = startState + alphabet_->getMin();
    int chrEnd = endState + alphabet_->getMin();
    double sumOfRates = 0; //for normalization of weights

    if (stateJumpTypeProb_.find(jumpStates) == stateJumpTypeProb_.end()){
        //gain
        if (chrStart + 1 == chrEnd){
            if (!(model->isIgnoredGain())){
                if (chrStart > 3){
                    stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::GAIN_T] = 1;
                    return;
                }else{
                    double gainRate = model->getGain()->getRate(chrStart);
                    stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::GAIN_T] = gainRate;
                    sumOfRates += gainRate;
                    legalMove = true;
                }

            }
        }
        //loss
        if (chrStart - 1 == chrEnd){
            if (!(model->isIgnoredLoss())){
                stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::LOSS_T] = 1;
                return;
            }
        }
        //baseNumber transitions
        int baseNumber = model->getBaseNumber();
        if (baseNumber != IgnoreParam){
            if (chrEnd > chrStart){
                if (((chrEnd - chrStart) % baseNumber == 0) && ((chrEnd - chrStart) <= (int)(model->getMaxChrRange()))){
                    legalMove = true;
                    stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::BASENUM_T] = model->getBaseNumR()->getRate(chrStart);
                    sumOfRates += model->getBaseNumR()->getRate(chrStart);                              
                }
            }        
        }
        //duplication
        if (chrEnd == 2 * chrStart){
            if (!(model->isIgnoredDupl())){
                legalMove = true;
                double duplRate = model->getDupl()->getRate(chrStart);
                stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::DUPL_T] = duplRate;
                sumOfRates += duplRate;
            }

        }
        //demi-duplication
        if (!(model->isIgnoredDemiDupl())){
            if (chrStart % 2 == 0){
                if (chrEnd == chrStart * 1.5){
                    legalMove = true;
                    stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::DEMIDUPL_T] = model->getDemiDupl()->getRate(chrStart);
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
                    stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::DEMIDUPL_T] = demiDupRate;
                    sumOfRates += demiDupRate;
                }
            }

        }

        if (chrEnd  == (int)(alphabet_->getMax())){
            legalMove = true;
            double toMaxRate = model->Qij(startState, endState)-sumOfRates;
            stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::MAXCHR_T] = toMaxRate;
            sumOfRates += toMaxRate;
        }
        //if nothing fits
        if(!legalMove){
            throw Exception ("ERROR: ComputeChromosomeTransitionsExp::updateMapOfJumps(): Illegal transition!");
            return;
        }
        // normalize according to weights
        map <int, double>::iterator it = stateJumpTypeProb_[jumpStates].begin();
        while (it != stateJumpTypeProb_[jumpStates].end()){
            int typeOfJump = it->first;
            double rate = stateJumpTypeProb_[jumpStates][typeOfJump];
            stateJumpTypeProb_[jumpStates][typeOfJump] = rate / sumOfRates;
            it ++;
        }
            
    }


}
// /********************************************************************************/
void ComputeChromosomeTransitionsExp::updateExpectationsPerBranch(uint nodeId, pair<int, int> ancestralTerminals, pair<int, int> jumpStates){
    map <int, double>::iterator it = stateJumpTypeProb_[jumpStates].begin();
    while (it != stateJumpTypeProb_[jumpStates].end()){
        int typeOfJump = it->first;
        double prob = stateJumpTypeProb_[jumpStates][typeOfJump];
        (branchTransitionsExp_[nodeId][ancestralTerminals].second)[typeOfJump] += prob;
        it ++;
    }
    
}

// /********************************************************************************/
int ComputeChromosomeTransitionsExp::getRandomState(int currentState, std::shared_ptr<const ChromosomeSubstitutionModel> model){
    double prob = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
    double cumulativeProb = 0;
    int nextState = currentState;
    for (size_t i = 0; i < alphabet_->getSize(); i++){
        double transitionRate;
        ((int)i == currentState) ? (transitionRate = 0) : (transitionRate = model->Qij(currentState,(int)i));
        cumulativeProb += transitionRate/(-1 * (model->Qij(currentState, currentState))); //jumpProbs_[currentState][(int)i];
        if (prob < cumulativeProb){
            nextState = (int)i;
            return nextState;
        }
    }
    throw Exception("ERROR: ComputeChangesExpectations::getRandomState(): could not sample new state.");
    return 1;
}
// /*************************************************************************************/
double ComputeChromosomeTransitionsExp::getCumulativeProbability(uint nodeId, vector <pair<int, int>>* terminalsToAccount){
    double cumulativeProb  = 0;
    std::map<std::pair<int, int>, std::pair<int, Vdouble>> :: iterator iterTerminalStates = branchTransitionsExp_[nodeId].begin();
    while(iterTerminalStates != branchTransitionsExp_[nodeId].end()){
        std::pair <int, int> currentPairOfAncestralTerminals = iterTerminalStates->first;
        // adding the accounted for terminals
        if(terminalsToAccount){
            terminalsToAccount->push_back(currentPairOfAncestralTerminals);

        }
        
        cumulativeProb += jointProbabilitiesFatherSon_[nodeId][0][currentPairOfAncestralTerminals.second][currentPairOfAncestralTerminals.first];
        iterTerminalStates ++;
    }
    return cumulativeProb;

}
// /*************************************************************************************/
double ComputeChromosomeTransitionsExp::isNeededHeuristics(uint nodeId, map <uint, vector<pair<int,int>>>* unAccountedNodesAndTerminals){
    vector <pair<int, int>> terminalsToAccount;
    double cumulativeProb = getCumulativeProbability(nodeId, &terminalsToAccount);
    if (cumulativeProb < THRESHOLD_HEURISTIC){
        (*unAccountedNodesAndTerminals)[nodeId] = terminalsToAccount;       
        

    }else{
        if (unAccountedNodesAndTerminals->find(nodeId) != unAccountedNodesAndTerminals->end()){
            unAccountedNodesAndTerminals->erase(nodeId);
        }

    }    
    return cumulativeProb;

}
// /*************************************************************************************/
void ComputeChromosomeTransitionsExp::updateNumNonAccountedBranches(map <uint, vector<pair<int,int>>>* unAccountedNodesAndTerminals, int iteration, size_t modelIndex, const string filePath){
    vector <double> cumulativeProbs;
    vector<Branch> branches;
    for (size_t n = 0; n < branchOrder_[modelIndex].size(); n++){
        Branch branchEdge(branchOrder_[modelIndex][n].first, branchOrder_[modelIndex][n].second);
        branches.push_back(branchEdge);
    }
   
    for (size_t i = 0; i < branches.size(); i++){
        uint nodeId = tree_->getSon(branches[i].first);
        double thresholdHeuristics = isNeededHeuristics(nodeId, unAccountedNodesAndTerminals);
        if (thresholdHeuristics >= THRESHOLD_HEURISTIC){
            Branch branchToRemove = branches[i];
            auto removed = remove_if(branchOrder_[modelIndex].begin(), branchOrder_[modelIndex].end(), [branchToRemove](const Branch& x) { return x.first == branchToRemove.first;});
            branchOrder_[modelIndex].erase(removed, branchOrder_[modelIndex].end());


        }else{
            if (iteration == 0){
                cumulativeProbs.push_back(thresholdHeuristics);
            }

        }
        

    }
    if (iteration == 0){
        
        if (filePath == "none"){
            return;
        }
        ofstream outFile;
        outFile.open(filePath);
        outFile << "Non Accounted for transitions:"<<endl;
        for (size_t k = 0; k < cumulativeProbs.size(); k++){
            string nodeName;
            uint nodeId = tree_->getSon(branchOrder_[modelIndex][k].first);
            //int nodeId = branchOrder_[k].getId();
            if (tree_->isLeaf(nodeId)){
                nodeName = (tree_->getNode(nodeId))->getName();

            }else{
                nodeName = "N" + std::to_string(nodeId);
            }
            outFile << nodeName <<": "<<cumulativeProbs[k]<<endl;

        }
        outFile.close();


    }


}
// /*************************************************************************************/
void ComputeChromosomeTransitionsExp::runHeuristics(const string FilePath){
    map <uint, vector<pair<int, int>>> nonAccountedForBranchesFromFirstRun;
    //map <int, double> branchMultiplier;
    map <uint, vector<pair<int,int>>> unAccountedNodesAndTerminals;
    map <int, double> ratesPerState;
    for (size_t m = 0; m < branchOrder_.size(); m ++){
        for (int i = 0; i < MAX_ITER_HEURISTICS; i++){
            updateNumNonAccountedBranches(&unAccountedNodesAndTerminals, i, m, FilePath);
            if (i == 0){
                nonAccountedForBranchesFromFirstRun = unAccountedNodesAndTerminals;
            }
            if (branchOrder_[m].size() == 0){
                break;
            }
            for (int k = 0; k < (int)(alphabet_->getSize()); k++){
                updateBranchLengths(k, i, m, &ratesPerState);
                for (int j = 0; j < MAX_SIM_HEURISTICS; j++){
                    runIteration(k, m, &unAccountedNodesAndTerminals);
                }
            }

        }

    }

    getPosteriorAndExpForNonAccountedFor(nonAccountedForBranchesFromFirstRun); 
    computeExpPerTypeHeuristics(nonAccountedForBranchesFromFirstRun);
}
// /*************************************************************************************/
void ComputeChromosomeTransitionsExp::updateBranchLengths(int initState, int iteration, size_t modelIndex, map <int, double>* ratesPerState){
    auto model = std::dynamic_pointer_cast<const ChromosomeSubstitutionModel>(model_->getModel(modelIndex+1));
    double sumOfRates = 0;
    if (iteration == 0){
        vector<double> rates;
        //dupl
        (model->isIgnoredDupl()) ? (rates.push_back(IgnoreParam)) : (rates.push_back(model->getDupl()->getRate(initState)));
        //demi-dupl
        (model->isIgnoredDemiDupl()) ? (rates.push_back(IgnoreParam)) : (rates.push_back(model->getDemiDupl()->getRate(initState)));
        //base number
        double baseNumRate = 0;
        int baseNumber = model->getBaseNumber();
        // we would like to take into consideration the overall rate for a base number transition from the current state
        if (baseNumber != IgnoreParam){
            for (int i = initState + baseNumber; i < model->getMax()-model->getMin(); i+=baseNumber){
                if ((i - initState) > (int)(model->getMaxChrRange())){
                    break;
                }
                baseNumRate += model->getBaseNumR()->getRate(initState);

            }
        }
        (baseNumber == IgnoreParam) ? (rates.push_back(IgnoreParam)) : (rates.push_back(baseNumRate));
        //gain
        (model->isIgnoredGain()) ? (rates.push_back(IgnoreParam)) : (rates.push_back(model->getGain()->getRate(initState)));
        //loss
        rates.push_back(model->getLoss()->getRate(initState));
        for (size_t i = 0; i < rates.size(); i++){
            if (rates[i] == IgnoreParam){
                continue;
            }
            sumOfRates += rates[i];

        }
        (*ratesPerState)[initState] = sumOfRates;

    }


    //multiplying by factor
    for (size_t k = 0; k < branchOrder_[modelIndex].size(); k++){
        if (iteration == 0){
            // the branch is too small, and the expected number of transitions is less than 1
            // Therefore, the branch length is scaled to be larger, so that the expected number of changes
            // will be at least 1
            if ((*ratesPerState)[initState] * branchOrder_[modelIndex][k].second.getLength() < 1){
                branchOrder_[modelIndex][k].second.setLength(1/(*ratesPerState)[initState]);
            }
        }else{
            // in each iteration, the extension factor is increased to increase the chance for transitions to occur
            if (tree_->getEdge(branchOrder_[modelIndex][k].first)->getLength() * (*ratesPerState)[initState] >= 1){
                branchOrder_[modelIndex][k].second.setLength(tree_->getEdge(branchOrder_[modelIndex][k].first)->getLength() * BRANCH_MULTIPLIER_FACTOR * iteration);

            }else{
                branchOrder_[modelIndex][k].second.setLength((1/(*ratesPerState)[initState]) * (BRANCH_MULTIPLIER_FACTOR * iteration));
            }
            
        }
         

    }
    sort(branchOrder_[modelIndex].begin(), branchOrder_[modelIndex].end(), compareBranches);

}

// /*************************************************************************************/
void ComputeChromosomeTransitionsExp::getPosteriorAndExpForNonAccountedFor(map <uint, vector<pair<int, int>>>& nonAccountedForBranchesFromFirstRun){
    map <uint, vector<pair<int, int>>>::iterator it = nonAccountedForBranchesFromFirstRun.begin();
    while (it != nonAccountedForBranchesFromFirstRun.end()){
        uint nodeId = it->first;
        auto itAllTerminals = branchTransitionsExp_[nodeId].begin();
        while(itAllTerminals != branchTransitionsExp_[nodeId].end()){
            int ancestralTerminal = itAllTerminals->first.first;
            int sonTerminal = itAllTerminals->first.second;
            std::pair<int, int> terminals(ancestralTerminal, sonTerminal);
            auto accountedForTerminals = std::find(nonAccountedForBranchesFromFirstRun[nodeId].begin(), nonAccountedForBranchesFromFirstRun[nodeId].end(), terminals);
            if (accountedForTerminals != nonAccountedForBranchesFromFirstRun[nodeId].end()){
                itAllTerminals ++;
                continue;

            }else{
                for (size_t j = 0; j < ChromosomeSubstitutionModel::NUMTYPES; j++){
                    (branchTransitionsExp_[nodeId][terminals].second)[j] /= branchTransitionsExp_[nodeId][terminals].first;
                }
                itAllTerminals ++; 

            }
           
        }

        it++;
    }

}

// /*************************************************************************************/
void ComputeChromosomeTransitionsExp::computeExpectationAndPosterior(){
    std::map <uint, std::map<std::pair<int, int>, pair<int, Vdouble>>>::iterator it = branchTransitionsExp_.begin();
    while (it != branchTransitionsExp_.end()){
        uint nodeId = it->first;
        std::map<std::pair<int, int>, pair<int, Vdouble>> :: iterator iterTerminalStates = branchTransitionsExp_[nodeId].begin();
        while(iterTerminalStates != branchTransitionsExp_[nodeId].end()){
            std::pair <int, int> currentPairOfAncestralTerminals = iterTerminalStates->first;
            int countForPairOfTerminals = branchTransitionsExp_[nodeId][currentPairOfAncestralTerminals].first;
            if (countForPairOfTerminals == 0){
                iterTerminalStates ++;
                continue;
            }
            for (size_t i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
                (branchTransitionsExp_[nodeId][currentPairOfAncestralTerminals].second)[i] /= branchTransitionsExp_[nodeId][currentPairOfAncestralTerminals].first;
            }

            iterTerminalStates ++;
        }        
        it ++;
    }
}
// /***************************************************************************************/
double ComputeChromosomeTransitionsExp::getExpectation(uint nodeId, int startAncestral, int endAncestral, int typeOfChange){
    std::pair <int, int> ancestralTerminals;
    ancestralTerminals.first = startAncestral;
    ancestralTerminals.second = endAncestral;
    return (branchTransitionsExp_[nodeId][ancestralTerminals].second)[typeOfChange];

}




/************************************************************************************************/
std::map<int, double> ComputeChromosomeTransitionsExp::getExpectationsPerType(const NonHomogeneousSubstitutionProcess* NonHomoProcess, PhyloTree &tree, std::map<uint, std::map<pair<size_t, size_t>, double>> &expectationsPerNode){
    std::map<int, double> expectationsPerType;
    std::map<uint, size_t> modelsForBranch = StochasticMappingUtils::getModelForEachBranch(tree, *NonHomoProcess); //son end of the branch and its corresponding model (father's model)
    auto it = expectationsPerNode.begin();
    while (it != expectationsPerNode.end()){
        auto nodeId = it->first;
        auto model = NonHomoProcess->getModel(modelsForBranch[nodeId]);
        auto chrModel = std::dynamic_pointer_cast<const ChromosomeSubstitutionModel>(model);
        if (expectationsPerNode.find(nodeId) == expectationsPerNode.end()){
            it ++;
            continue;
        }
        auto &transitionsPerNode = expectationsPerNode[nodeId]; // don't want to create a local copy of this element, just to use a reference
        auto transitionsTypesPerNode = StochasticMappingUtils::getTypeForEachTransitionPerNode(chrModel, transitionsPerNode, nodeId);
        auto typesIt = transitionsTypesPerNode.begin();
        while(typesIt != transitionsTypesPerNode.end()){
            if (transitionsTypesPerNode.find(typesIt->first) == transitionsTypesPerNode.end()){
                expectationsPerType[typesIt->first] = transitionsTypesPerNode[typesIt->first];
            }else{
                expectationsPerType[typesIt->first] += transitionsTypesPerNode[typesIt->first];
            }
            typesIt ++;
        }
       
        it ++;
    }
    return expectationsPerType;
}

/************************************************************************************/
