#include "ChromosomeNumberOptimizer.h"
using namespace bpp;


/**********************************************************************************/
void ChromosomeNumberOptimizer::fillVectorOfLikelihoods(SingleProcessPhyloLikelihood* lik, uint numOfIterationsFirstCycle,  size_t currPoint, uint reqNumOfPoints, vector <uint> baseNumCandidates, std::map<int, vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>& fixedParams, vector<SingleProcessPhyloLikelihood*> &vectorOfLiklihoods, string* text, std::map<uint, uint> &baseNumberUpperBounds, omp_lock_t* mutex){
    printLikParameters(lik, 0, text);
    if (numOfIterationsFirstCycle > 0){
        optimizeModelParameters(lik, tolerance_, numOfIterationsFirstCycle, baseNumCandidates, sharedParams, &fixedParams, text, baseNumberUpperBounds);
        printLikParameters(lik, 1, text);

    }
    if(mutex){
        omp_set_lock(mutex);
        //std::cout << "Within critiacal section: fillVectorOfLikelihoods()" << std::endl;

    }
    if (currPoint < reqNumOfPoints){
        vectorOfLiklihoods.push_back(lik);
            
    }else{
        SingleProcessPhyloLikelihood* likToDel;
        size_t index;

        for (size_t i = 0; i < reqNumOfPoints; i++){ 
            if (i == 0){
                likToDel = vectorOfLiklihoods[0];
                index = 0;
            }else{
                if (likToDel->getValue() < vectorOfLiklihoods[i]->getValue()){
                    likToDel = vectorOfLiklihoods[i];
                    index = i;
                }
            }
        }
        if (lik->getValue() < likToDel->getValue()){
            vectorOfLiklihoods[index] = lik;
            

        }else{
            likToDel = lik;

        }

        deleteLikObject(likToDel);
                      

    }
    if (mutex){
        //std::cout << "Out of critiacal section: fillVectorOfLikelihoods()" << std::endl;
        omp_unset_lock(mutex);
            

    }      

}
/**********************************************************************************/
void ChromosomeNumberOptimizer::initLikelihoods(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds, uint numOfModels, std::map<int, vector<std::pair<uint, int>>>* sharedParams){
    size_t index = min((int)numOfPoints_.size()-1, 1);
    vectorOfLikelohoods_.reserve(numOfPoints_[index]);
    vector <unsigned int> baseNumCandidates;
    getBaseNumCandidates(baseNumCandidates, baseNumberUpperBound_);

    // If base number is one of the parameters
    SingleProcessPhyloLikelihood* lik;
    cout <<"##################################" << endl;
    cout << "*********  cycle 0  **************"<<endl;  
    for (size_t n = 0; n < numOfPoints; n++){
        std::cout << "Starting cycle with Point #" << n <<"...."<<endl;
        if (n == 0){
            lik = setHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels, sharedParams);
        }else{
            if (n == 1){
                lik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels, parsimonyBound * (double)n, fixedParams, sharedParams);

            }else{
                lik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels, parsimonyBound * (1+(0.1*(double)n)), fixedParams, sharedParams);
            }
            
    
        }
        ifNanTryToResampleLikObject(&lik, tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels, parsimonyBound, numOfPoints, fixedParams, sharedParams);
        // int numOfTrials
        
        // while((std::isnan(lik->getValue())) && (ChromEvolOptions::maxNumOfTrials_)){

        // }
        fillVectorOfLikelihoods(lik, numOfIterations_[0],  n, numOfPoints_[index], baseNumCandidates, sharedParams, fixedParams, vectorOfLikelohoods_, 0, baseNumberUpperBound_);
   
    }

    sort(vectorOfLikelohoods_.begin(), vectorOfLikelohoods_.end(), compareLikValues);
    printLikelihoodVectorValues(vectorOfLikelohoods_, 0, 0);

}

// /****************************************************************************/
vector <double> ChromosomeNumberOptimizer::setFixedRootFrequencies(const std::string &path, std::shared_ptr<ChromosomeSubstitutionModel> chrModel){
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

/*******************************************************************************/
void ChromosomeNumberOptimizer::getBaseNumCandidates(vector <unsigned int> &baseNumCandidates, std::map<uint, uint> &baseNumberUpperBounds) const{
    if ((baseNumOptimizationMethod_ != "Brent") && (optimizeBaseNumber_)){
        uint maxBaseNumCandidate = getMaxBaseNumAmongModels(baseNumberUpperBounds);
        fillVectorOfBaseNumCandidates(baseNumCandidates, lowerBoundBaseNumber, maxBaseNumCandidate);

    }

}

// /****************************************************************************/

void ChromosomeNumberOptimizer::optimizeMultiProcessModel(std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, vector<unsigned int> &numOfPoints, vector<unsigned int> &numOfIterations, std::map<uint, uint> &baseNumberUpperBounds, vector<SingleProcessPhyloLikelihood*>* perCandidateLik, string* text, omp_lock_t* mutex)
{
    if (mutex){
        omp_set_lock(mutex);
        //std::cout << "Within critical section: optimizeMultiProcessModel(): base number fill" << std::endl;

    }
    vector <unsigned int> baseNumCandidates;
    getBaseNumCandidates(baseNumCandidates, baseNumberUpperBounds);
    if (mutex){
        //std::cout << "Out of critical section: optimizeMultiProcessModel():base number fill" << std::endl;
        omp_unset_lock(mutex);
    }
    //Go over each cycle
    // note: I start from the second point, because the first iteration has been already done.
    for (size_t i = 1; i < numOfIterations.size(); i++){
        for (size_t j = 0; j < numOfPoints[i]; j++){
            if (!perCandidateLik){
                printLikParameters(vectorOfLikelohoods_[j], 0, text);
            }
            //If the number of optimization iterations is larger than zero, optimize the number of times as specified
            if (numOfIterations[i] > 0){
                if (perCandidateLik){
                    // add here baseNum candidates
                    optimizeModelParameters((*perCandidateLik)[j], tolerance_, numOfIterations[i], baseNumCandidates, sharedParams, fixedParams, text, baseNumberUpperBounds);
                }else{
                    optimizeModelParameters(vectorOfLikelohoods_[j], tolerance_, numOfIterations[i], baseNumCandidates, sharedParams, fixedParams, text, baseNumberUpperBounds);
                }
            }
            if ((i < numOfIterations.size()-1) && (j >= numOfPoints[i+1])){
                if (mutex){
                    omp_set_lock(mutex);
                    //std::cout << "Within critical section: optimizeMultiProcessModel()" << std::endl;

                }
                
                SingleProcessPhyloLikelihood* currLik;
                SingleProcessPhyloLikelihood* likToDel;
                (perCandidateLik) ? (currLik = (*perCandidateLik)[j]) : (currLik = vectorOfLikelohoods_[j]);
                size_t index = 0;
                for (size_t k = 0; k < numOfPoints[i+1]; k++){
                    if (k == 0){
                        if (perCandidateLik){
                            likToDel = (*perCandidateLik)[k];

                        }else{
                            likToDel = vectorOfLikelohoods_[k];
                        }
                        
                    }else{
                        if (perCandidateLik){
                            if (likToDel->getValue() < (*perCandidateLik)[k]->getValue()){
                                likToDel = (*perCandidateLik)[k];
                                index = k;
                            }
                        }else{
                            if (likToDel->getValue() < vectorOfLikelohoods_[k]->getValue()){
                                likToDel = vectorOfLikelohoods_[k];
                                index = k;
                            }

                        }
                    }

                }
                if (currLik->getValue() > likToDel->getValue()){
                    likToDel = currLik;
                }else{
                    (perCandidateLik) ? ((*perCandidateLik)[index] = currLik) : (vectorOfLikelohoods_[index] = currLik);
                }
                
                deleteLikObject(likToDel);
                (perCandidateLik) ? ((*perCandidateLik)[j] = 0) : (vectorOfLikelohoods_[j] = 0);
                if (mutex){
                    //std::cout << "Out of critical section: optimizeMultiProcessModel()" << std::endl;
                    omp_unset_lock(mutex);
                }
            }     
        }
        if (i < numOfIterations.size()-1){
            if (mutex){
                omp_set_lock(mutex);
                //std::cout << "Within critical section2 : optimizeMultiProcessModel()" << std::endl;
            }
            if (perCandidateLik){
                clearVectorOfLikelihoods(numOfPoints[i+1], *perCandidateLik);

            }else{
                clearVectorOfLikelihoods(numOfPoints[i+1]);
            }
            if (mutex){
                //std::cout << "Out of critical section2: optimizeMultiProcessModel()" << std::endl;
                omp_unset_lock(mutex);
            }
                      
        }
        if (mutex){
            omp_set_lock(mutex);
            //std::cout << "Within critical section3: optimizeMultiProcessModel()" << std::endl;
        }

        if (perCandidateLik){
            sort((*perCandidateLik).begin(), (*perCandidateLik).end(), compareLikValues);
            printLikelihoodVectorValues(*perCandidateLik, text, i);
        }else{
            sort(vectorOfLikelohoods_.begin(), vectorOfLikelohoods_.end(), compareLikValues);
            printLikelihoodVectorValues(vectorOfLikelohoods_, 0, i);

        }
        if (mutex){
            //std::cout << "Out of critical section3: optimizeMultiProcessModel()" << std::endl;
            omp_unset_lock(mutex);
        }
    }


}
/********************************************************************************/
void ChromosomeNumberOptimizer::deleteLikObject(SingleProcessPhyloLikelihood* lik_to_del){
    auto sequenceData = lik_to_del->getData();
    auto process = &(lik_to_del->getSubstitutionProcess());
    //auto tree = &(lik_to_del->getTree());
    auto context = &(lik_to_del->getContext());
    delete process;
    delete sequenceData;
    //delete tree;
    delete context;
    delete lik_to_del;
    

}

// /********************************************************************************/
void ChromosomeNumberOptimizer::clearVectorOfLikelihoods(size_t new_size){
    while(vectorOfLikelohoods_.size() > new_size){
        //deleteTreeLikAssociatedAttributes(vectorOfLikelohoods_[vectorOfLikelohoods_.size()-1]);
        SingleProcessPhyloLikelihood* lik_to_del = vectorOfLikelohoods_.back(); 
        vectorOfLikelohoods_.pop_back();
        if (lik_to_del){
            deleteLikObject(lik_to_del);

        }
    }
}
// /********************************************************************************/
void ChromosomeNumberOptimizer::clearVectorOfLikelihoods(size_t new_size, std::vector<SingleProcessPhyloLikelihood*> &likelihoodsVec){
    while(likelihoodsVec.size() > new_size){
        //deleteTreeLikAssociatedAttributes(vectorOfLikelohoods_[vectorOfLikelohoods_.size()-1]);
        SingleProcessPhyloLikelihood* lik_to_del = likelihoodsVec.back(); 
        likelihoodsVec.pop_back();
        if (lik_to_del){
            deleteLikObject(lik_to_del);

        }
        

    }
}

// /***********************************************************************************/
bool ChromosomeNumberOptimizer::compareLikValues(SingleProcessPhyloLikelihood* lik1, SingleProcessPhyloLikelihood* lik2){
    return (lik1->getValue() < lik2->getValue());
}
// /***********************************************************************************/
void ChromosomeNumberOptimizer::printLikParameters(SingleProcessPhyloLikelihood* lik, unsigned int optimized, string* textToPrint, const string filePath) const{
    ofstream outFile;
    if (filePath != "none"){
        outFile.open(filePath);
    }
    string text;
    if (optimized == 0){
        text = "Initial likelihood is : "+ std::to_string(lik->getValue())+"\n";
    }else{
        text = "Optimized likelihood is : "+ std::to_string(lik->getValue()) +"\n";
        if (filePath != "none"){
            outFile << "Final optimized likelihood is: "<< lik->getValue() << endl;
        }
    }
    text +=  "Parameters are:\n";
    if (filePath != "none"){
        outFile << "Optimized parameters are:"<<endl;
    }
    ParameterList substitutionModelParams = lik->getSubstitutionModelParameters();
    std::vector<std::string> paramsNames = substitutionModelParams.getParameterNames();
    for (int i = 0; i < (int)(paramsNames.size()); i++){
        if (paramsNames[i].find("Chromosome.baseNum_") != std::string::npos){
            text += paramsNames[i]+ " = "+ std::to_string((int)(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()))+"\n";
            if (filePath != "none"){
                outFile <<  paramsNames[i] << " = "<< (int)(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()) <<endl;
            }
        }else{
            text += (paramsNames[i] + " = "+ std::to_string(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()) + "\n");
            if (filePath != "none"){
                outFile << paramsNames[i] << " = "<< lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue() <<endl;
            }
        }
        
    }
    if (filePath != "none"){
        outFile.close();
    }
    text +=  "***\n";
    if (textToPrint){
        *textToPrint += text;
    }else{
        std::cout << text;
    }

}
// /*************************************************************************************/
void ChromosomeNumberOptimizer::printLikelihoodVectorValues(std::vector <SingleProcessPhyloLikelihood*> lik_vec, string* text, size_t cycle) const{
    string textToPrint = "The likelihoods at the end of cycle " + std::to_string(cycle)+" are :\n";
    for (size_t i = 0; i < lik_vec.size(); i++){
        textToPrint += std::to_string(lik_vec[i]->getValue())+"\n";
    }
    if (text){
        *text += textToPrint;
    }else{
        std::cout << textToPrint;
    }
}
/********************************************************************************/
std::vector<double> ChromosomeNumberOptimizer::getRootFrequencies(SingleProcessPhyloLikelihood* lik) const{
    std::vector<double> rootFreqs;
    ValueRef <Eigen::RowVectorXd> rootFreqVector = lik->getLikelihoodCalculationSingleProcess()->getRootFreqs();
    for (size_t s = 0; s < (size_t)rootFreqVector->getTargetValue().size(); s ++){
        rootFreqs.push_back(rootFreqVector.get()->getTargetValue()[s]);

    }
    return rootFreqs;

}
// /******************************************************************************/
void ChromosomeNumberOptimizer::printRootFrequencies(SingleProcessPhyloLikelihood* lik, ofstream &outFile) const{
    ValueRef <Eigen::RowVectorXd> rootFreqVector = lik->getLikelihoodCalculationSingleProcess()->getRootFreqs();
    for (size_t s = 0; s < (size_t)rootFreqVector->getTargetValue().size(); s ++){
        cout << "F[" << s + alphabet_->getMin() << "] = " << rootFreqVector.get()->getTargetValue()[s] << endl;
        if (s == (size_t)(rootFreqVector->getTargetValue().size()) -1){
            outFile << "F[" << s + alphabet_->getMin() << "] = " << rootFreqVector.get()->getTargetValue()[s] << endl;
        }else{
            outFile << "F[" << s + alphabet_->getMin() << "] = " << rootFreqVector.get()->getTargetValue()[s] << "\t";

        }         
        
    }
}
/**************************************************************************************/
uint ChromosomeNumberOptimizer::getMaxBaseNumAmongModels(std::map<uint, uint> baseNumberUpperBound) const{
    auto it = baseNumberUpperBound.begin();
    uint maxBaseNumBound = 0;
    while(it != baseNumberUpperBound.end()){
        if (baseNumberUpperBound[it->first] > maxBaseNumBound){
            maxBaseNumBound = baseNumberUpperBound[it->first];
        }
        it ++;
    }
    return maxBaseNumBound;
}

// /***********************************************************************************/
void ChromosomeNumberOptimizer::fillVectorOfBaseNumCandidates(std::vector <unsigned int> &baseNumCandidates, unsigned int lowerBound, unsigned int upperBound) const{
    if (baseNumOptimizationMethod_ == "Ranges"){
        getAllPossibleChrRanges(baseNumCandidates);

    }
    else if ((baseNumOptimizationMethod_ == "Sequential") || (baseNumCandidates.size() == 0)){

        for (unsigned int chr = (unsigned int)lowerBound; chr <= upperBound; chr++){
            baseNumCandidates.push_back(chr);
        }

    }

}
// /***************************************************************************************/
void ChromosomeNumberOptimizer::getAllPossibleChrRanges(std::vector <unsigned int> &baseNumCandidates) const{
    size_t numOfSequences = vsc_->getNumberOfSequences();
    unsigned int minRange = 0;
    vector <string> sequenceNames = vsc_->getSequenceNames();
    for (size_t i = 0; i < numOfSequences; i++){
        if (i == numOfSequences-1){
            continue;
        }
        BasicSequence seq1 = vsc_->getSequence(sequenceNames[i]);
        int chrNum1 = seq1.getValue(0);
        if (chrNum1 == -1){
            continue;
        }
        for (size_t j = i + 1; j < numOfSequences; j++){
            BasicSequence seq2 = vsc_->getSequence(sequenceNames[j]);
            int chrNum2 = seq2.getValue(0);
            if (chrNum2 == -1){
                continue;
            }
            unsigned int chrRange = (unsigned int)(std::abs(chrNum1 - chrNum2));
            if (chrRange < lowerBoundBaseNumber){
                continue;
            }

            if (!std::count(baseNumCandidates.begin(), baseNumCandidates.end(), chrRange)){
                if (minRange == 0){
                    minRange = chrRange;
                }else{
                    if (chrRange < minRange){
                        minRange = chrRange;
                    }
                }
                baseNumCandidates.push_back(chrRange);

            }

        }
    }
    if (minRange > lowerBoundBaseNumber){
        for (unsigned int i = lowerBoundBaseNumber; i < minRange; i++){
            baseNumCandidates.push_back(i);
        }

    }

}

// /**********************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeModelParameters(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds){
    unsigned int numOfEvaluations = 0;
 
    if (typeOfOptimizer_ == "Brent"){
        numOfEvaluations += optimizeModelParametersOneDimension(tl, tol, maxNumOfIterations, baseNumCandidates, sharedParams, fixedParams, textToPrint, baseNumberUpperBounds);
    }else if (typeOfOptimizer_ == "gradient"){
        checkLegalUseOfGradientOptimization();
        numOfEvaluations += optimizeMultiDimensions(tl, tol, maxNumOfIterations, sharedParams, fixedParams, textToPrint, baseNumberUpperBounds);

    }else{
        checkLegalUseOfGradientOptimization();
        numOfEvaluations += useMixedOptimizers(tl, tol, maxNumOfIterations, baseNumCandidates, sharedParams, fixedParams, textToPrint, baseNumberUpperBounds);
    }
        
    return numOfEvaluations;
    
}
// /***************************************************************************************/
void ChromosomeNumberOptimizer::checkLegalUseOfGradientOptimization(){
    // Only base num is not a composite parameter
    size_t startForComposite = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();
    for (size_t k = startForComposite; k < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; k++){
        if (ChromEvolOptions::rateChangeType_[k-startForComposite] == ChromosomeNumberDependencyFunction::CONSTANT){
            continue;
        }else if (ChromEvolOptions::rateChangeType_[k-startForComposite] == ChromosomeNumberDependencyFunction::IGNORE){
            continue;
        }else if (ChromEvolOptions::rateChangeType_[k-startForComposite] == ChromosomeNumberDependencyFunction::EXP){
            continue;
        }else{
            throw Exception ("ChromosomeNumberOptimizer::checkLegalUseOfGradientOptimization: Cannot use a gradient descent optimization! Use only Brent!!!");

        }
    }
}

// /****************************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeMultiDimensions(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, bool mixed, unsigned int currentIterNum){
    DerivableSecondOrder* f = tl;
    ParameterList tmp = tl->getSubstitutionModelParameters();
    unique_ptr<AbstractNumericalDerivative> fnum;
    fnum.reset(new TwoPointsNumericalDerivative(f));
    fnum->setInterval(0.0000001);
    ConjugateGradientMultiDimensions* optimizer = new ConjugateGradientMultiDimensions(fnum.get());
    fnum->setParametersToDerivate(tmp.getParameterNames());
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->getStopCondition()->setTolerance(tol* 0.1);
    optimizer->setMaximumNumberOfEvaluations(1000);
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, tl);
    size_t startCompositeParams = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();

    string text = "";
    unsigned int numOfEvaluations = 0;
    double currentLikelihood = tl->getValue();
    double prevLikelihood;
    for (size_t i = 0; i < maxNumOfIterations; i++){
        if(mixed){
            text += ("Iteration #"+ std::to_string(currentIterNum)+"\n");

        }else{
            text += ("Iteration #"+ std::to_string(i)+"\n");
        }
        
        ParameterList paramsFull = tl->getSubstitutionModelParameters();
        std::vector <string> nonFixedparamsNames = getNonFixedParams(tl, paramsFull, fixedParams);
        ParameterList params = tl->getParameters().createSubList(nonFixedparamsNames);
        int rateParamType;
        double lowerBound;
        double upperBound;
        
        for (size_t j = 0; j < params.size(); j++){
            std::string nameOfParam = params[j].getName();
            rateParamType = paramNameAndType[nameOfParam].first;
            /////////////////////////////////////////////////////
            std::vector<string> paramsNames = typeWithParamNames[rateParamType][paramNameAndType[nameOfParam].second];
            auto it = std::find(paramsNames.begin(), paramsNames.end(), nameOfParam);
            if (it == paramsNames.end()){
                throw Exception("ChromosomeNumberOptimizer::optimizeModelParametersOneDimension(): index out of range!");
            }
            size_t index = it - paramsNames.begin();
            ////////////////////////////////////////////////////////
            if (rateParamType != ChromosomeSubstitutionModel::BASENUM){
                ChromosomeNumberDependencyFunction::FunctionType funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(ChromEvolOptions::rateChangeType_[rateParamType-startCompositeParams]);
                ChromosomeNumberDependencyFunction* functionOp = compositeParameter::setDependencyFunction(funcType);
                functionOp->setDomainsIfNeeded(alphabet_->getMin(), alphabet_->getMax());
                functionOp->updateBounds(params, paramsNames, index, &lowerBound, &upperBound, alphabet_->getMax());
                std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(nameOfParam).getConstraint());
                interval->setLowerBound(lowerBound, interval->strictLowerBound());
                functionOp->updateBounds(f, nameOfParam, lowerBound, upperBound);
                delete functionOp;

            }    

        }
        prevLikelihood = currentLikelihood;
        optimizer->init(params);
        currentLikelihood = optimizer->optimize();
        *textToPrint += text;
        if (!textToPrint){
            std::cout << text;
        }
        printLikParameters(tl, 1, textToPrint);
        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        
        
    }
    
    numOfEvaluations += optimizer->getNumberOfEvaluations();
    if (!mixed){
        if (textToPrint){
            *textToPrint += "...\n";
        }else{
            std::cout <<"..."<<endl;

        }
    }
    //std::cout << "The final number of evaluations is: "<< numOfEvaluations << endl;
    delete optimizer;
    return numOfEvaluations;

}
// /*******************************************************************************/

unsigned int ChromosomeNumberOptimizer::useMixedOptimizers(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, vector <unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds){
    std::vector<size_t> optimization = RandomTools::randMultinomial(maxNumOfIterations, probsForMixedOptimization_);
    unsigned int numOfEvaluations = 0;
    for (size_t i = 0; i < maxNumOfIterations; i++){
        double prevLikelihood = tl->getValue();
        if (optimization[i] == 0){
            string text = "Optimizing with Brent\n";
            printLog(textToPrint, text);
            numOfEvaluations += optimizeModelParametersOneDimension(tl, tol, 1, baseNumCandidates, sharedParams, fixedParams, textToPrint, baseNumberUpperBounds, true, (unsigned int)i);
        }else{
            string text = "Optimizing with Gradient Descent\n";
            printLog(textToPrint, text);
            numOfEvaluations += optimizeMultiDimensions(tl, tol, 1, sharedParams, fixedParams, textToPrint, baseNumberUpperBounds, true, (unsigned int)i);
        }
        double currentLikValue = tl->getValue();
        if (std::abs(prevLikelihood-currentLikValue) < tol){
            break;
        }


    }
    
    return numOfEvaluations;

}
// /*******************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates, std::map<int, vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberBounds, bool mixed, unsigned curentIterNum){

    // Initialize optimizer
    string text;
    DerivableSecondOrder* f = tl;
    BrentOneDimension* optimizer = new BrentOneDimension(f);
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);
    std::cout <<"max chromosome number: " << alphabet_->getMax() << endl;
    size_t startCompositeParams = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();

    // setting bracketing for Brent optimization
    if (BrentBracketing_ == 1){
        optimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);

    }else if (BrentBracketing_ == 2){
        optimizer->setBracketing(BrentOneDimension::BRACKET_SIMPLE);
    }else{
        optimizer->setBracketing(BrentOneDimension::BRACKET_OUTWARD);
    }
    // initializing the likelihood values
    double currentLikelihood = tl->getValue();
    double prevLikelihood;
    unsigned int numOfEvaluations = 0;
    // setting maps of parameter type and the corresponding parameters, and vice versa
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, tl);
    ParameterList params;
    // starting iterations of optimization
    for (size_t i = 0; i < maxNumOfIterations; i++){
        if (mixed){
            text = "Iteration #" + std::to_string(curentIterNum)+ "\n";
            

        }else{
            text = "Iteration #" +std::to_string(i) + "\n";
        }
        printLog(textToPrint, text);
        //ParameterList params = tl->getParameters();// = tl->getParameters();
        ParameterList substitutionModelParams = tl->getSubstitutionModelParameters();
        size_t nbParams = substitutionModelParams.size();
        prevLikelihood = currentLikelihood;
        
        for (size_t j = 0; j < nbParams; j ++){
            params = tl->getParameters();
            const string nameOfParam = substitutionModelParams[j].getName();
            text = "Previous value of "+ nameOfParam + " is: "+ std::to_string(params.getParameter(nameOfParam).getValue()) + "\n";
            printLog(textToPrint, text);
            int rateParamType = paramNameAndType[nameOfParam].first;
            
            if (std::count((*fixedParams)[paramNameAndType[nameOfParam].second].begin(), (*fixedParams)[paramNameAndType[nameOfParam].second].end(), rateParamType)){
                continue;
            }
            
            //int rateCompositeParamType;
            double lowerBound;
            double upperBound;
            // param names corresponding to the parameter type
            std::vector<string> paramsNames = typeWithParamNames[rateParamType][paramNameAndType[nameOfParam].second];
            Parameter param = params.getParameter(nameOfParam);


            auto it = std::find(paramsNames.begin(), paramsNames.end(), nameOfParam);
            if (it == paramsNames.end()){
                throw Exception("ChromosomeNumberOptimizer::optimizeModelParametersOneDimension(): index out of range!");
            }
            size_t index = it - paramsNames.begin();
            if (rateParamType != static_cast<int>(ChromosomeSubstitutionModel::BASENUM)){
                ChromosomeNumberDependencyFunction::FunctionType funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(ChromEvolOptions::rateChangeType_[rateParamType-startCompositeParams]);
                ChromosomeNumberDependencyFunction* functionOp = compositeParameter::setDependencyFunction(funcType);
                functionOp->setDomainsIfNeeded(alphabet_->getMin(), alphabet_->getMax());
                functionOp->updateBounds(params, paramsNames, index, &lowerBound, &upperBound, alphabet_->getMax());
                functionOp->updateBounds(f, nameOfParam, lowerBound, upperBound);
                delete functionOp;

                std::shared_ptr<IntervalConstraint> intervalFuncUpdated = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(nameOfParam).getConstraint());
                //double updated_lowerBound = intervalFuncUpdated->getLowerBound();
                //std::cout << "*** ***" << nameOfParam << ": Updated lower bound: " << updated_lowerBound << std::endl;


                std::shared_ptr<IntervalConstraint> intervalFuncUpdatedTL = dynamic_pointer_cast<IntervalConstraint>(tl->getParameter(nameOfParam).getConstraint());
                //double updated_lowerBoundTL = intervalFuncUpdatedTL->getLowerBound();
                //std::cout << "*** ***" << nameOfParam << ": Updated lower bound TL: " << updated_lowerBoundTL << std::endl;  

            }else{
                // baseNumber parameter
                if (baseNumOptimizationMethod_ != "Brent"){
                    if (!std::count((*fixedParams)[paramNameAndType[nameOfParam].second].begin(), (*fixedParams)[paramNameAndType[nameOfParam].second].end(), ChromosomeSubstitutionModel::BASENUM)){
                        optimizeBaseNum(tl, j, baseNumCandidates, &currentLikelihood, lowerBound, upperBound, nameOfParam, params, paramNameAndType[nameOfParam].second, baseNumberBounds);
                        text = "parameter value after optimization "+ std::to_string(tl->getLikelihoodCalculation()->getParameter(param.getName()).getValue())+ "\n";
                        printLog(textToPrint, text);
                        continue;
                    }
                }
            }
            text = "Parameter name is: "+ nameOfParam +"\n";
            printLog(textToPrint, text);         
                 
            if ((i == 1) & (maxNumOfIterations > 2)){
                optimizer->getStopCondition()->setTolerance(tol* 2);
            }else{
                optimizer->getStopCondition()->setTolerance(tol);
            }
            if (rateParamType != static_cast<int>(ChromosomeSubstitutionModel::BASENUM)){
                optimizer->setInitialInterval(lowerBound + 1e-10, upperBound);
            }else{
                optimizer->setInitialInterval(lowerBound, upperBound);
            }            
            optimizer->init(params.createSubList(param.getName()));
            currentLikelihood = optimizer->optimize();
            text = "parameter value after optimization "+ std::to_string(tl->getLikelihoodCalculation()->getParameter(param.getName()).getValue())+ "\n";
            printLog(textToPrint, text);
            text = "***\n";
            printLog(textToPrint, text);
        }
        printLikParameters(tl, 1, textToPrint);
        
        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        numOfEvaluations += optimizer->getNumberOfEvaluations();
       
    }
    if (!mixed){
        text =  "...\n";
        printLog(textToPrint, text);
    }
    delete optimizer;
    return numOfEvaluations;
}
/**********************************************************************************/
void ChromosomeNumberOptimizer::createMapOfSharedParameterNames(std::map<int, std::vector<std::pair<uint, int>>> &sharedParams, std::map<string, vector<std::pair<uint, int>>> &sharedParamsNames){
    auto it = sharedParams.begin();
    while (it != sharedParams.end()){
        auto sharedPerParamNum = sharedParams[it->first];
        uint model = sharedPerParamNum[0].first;
        int type = sharedPerParamNum[0].second;
        string paramBasicName = getStringParamName(type);
        string paramName;
        size_t numOfParams = getNumberOfParametersPerParamType(type, ChromEvolOptions::rateChangeType_);
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
/**********************************************************************************/
uint ChromosomeNumberOptimizer::getModelFromParamName(string name){
    std::regex modelPattern ("_([\\d]+)");
    std::smatch sm;
    std::regex_search(name, sm, modelPattern);
    std::string modelSuffix = sm[sm.size()-1];
    uint modelId = static_cast<uint>(stoi(modelSuffix));
    return modelId;

}
/**********************************************************************************/
int ChromosomeNumberOptimizer::getTypeOfParamFromParamName(string name){
    int type;
    std::map<std::string, int> typeGeneralName;
    updateWithTypeAndCorrespondingName(typeGeneralName);
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
// /*******************************************************************************/
void ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>>* paramNameAndType, SingleProcessPhyloLikelihood* tl, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams){
    ParameterList substitutionModelParams = tl->getSubstitutionModelParameters();
    std::vector<std::string> namesAllParams = substitutionModelParams.getParameterNames();
    std::map<string, vector<std::pair<uint, int>>> sharedParamsNames;
    if (sharedParams){
        createMapOfSharedParameterNames(*sharedParams, sharedParamsNames);
    }
    for (size_t i = 0; i < namesAllParams.size(); i++){
        uint modelId = getModelFromParamName(namesAllParams[i]);
        int type = getTypeOfParamFromParamName(namesAllParams[i]);
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
// /*******************************************************************************/
void ChromosomeNumberOptimizer::updateWithTypeAndCorrespondingName(std::map<std::string, int> &typeGeneralName){
    typeGeneralName["gain"] = static_cast<int>(ChromosomeSubstitutionModel::GAIN);
    typeGeneralName["loss"] = static_cast<int>(ChromosomeSubstitutionModel::LOSS);
    typeGeneralName["dupl"] = static_cast<int>(ChromosomeSubstitutionModel::DUPL);
    typeGeneralName["demi"] = static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL);
    typeGeneralName["baseNumR"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUMR);
    typeGeneralName["baseNum_"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUM);
    
}
/*******************************************************************************/
size_t ChromosomeNumberOptimizer::getNumberOfFixedParams(SingleProcessPhyloLikelihood* lik, std::map<uint, vector<int>> &fixedParams){
    auto substitutionModelParams = lik->getSubstitutionModelParameters();
    auto paramNames = substitutionModelParams.getParameterNames();
    size_t numOfFixedParams = 0;
    for (size_t i = 0; i < paramNames.size(); i++){
        auto paramName = paramNames[i];
        uint model = getModelFromParamName(paramName);
        int type = getTypeOfParamFromParamName(paramName);
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
/*******************************************************************************/
std::string ChromosomeNumberOptimizer::getFunctionName(int func){
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
        throw Exception("ChromEvolOptions::getFunctionName: parameter not found !!!");
    }
    return functionName;
}
/*******************************************************************************/
string ChromosomeNumberOptimizer::getStringParamName(int type){
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
        throw Exception("ChromosomeNumberOptimizer::getStringParamName(): No such parameter exists!");
    }
    return strName;
}
/*******************************************************************************/
int ChromosomeNumberOptimizer::getEnumOfParamName(std::string pattern){

    if (pattern == "gain"){
        return static_cast<int>(ChromosomeSubstitutionModel::GAIN);
    }else if (pattern == "loss"){
        return static_cast<int>(ChromosomeSubstitutionModel::LOSS);
    }else if (pattern == "dupl"){
        return static_cast<int>(ChromosomeSubstitutionModel::DUPL);
    }else if (pattern ==  "baseNumR"){
        return static_cast<int>(ChromosomeSubstitutionModel::BASENUMR);
    }else if (pattern ==  "baseNum_"){
        return static_cast<int>(ChromosomeSubstitutionModel::BASENUM);
    }else if (pattern ==  "demi"){
        return static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL);
    
    }else{
        throw Exception("ChromosomeNumberOptimizer::getEnumOfParamName(): No such parameter exists!");
    }
    return -1;

}

/*******************************************************************************/
std::map<uint, std::vector<string>> ChromosomeNumberOptimizer::getRelatedParameterNamesForEachModel(ParameterList &params, std::string pattern, uint numOfModels, std::map<int, vector<std::pair<uint, int>>>* mapSharedParams){
  std::map<uint, std::vector<string>> matchingParamsPerModel;
  for (uint i = 1; i <= numOfModels; i++){
    matchingParamsPerModel[i] = vector<string>();
  }
  std::vector<std::string> paramNames = params.getParameterNames();
  std::regex modelPattern ("_([\\d]+)");
  for (size_t i = 0; i < paramNames.size(); i++){
    std::string fullParamName = paramNames[i];
    if (fullParamName.find(pattern) != string::npos){
      // put to the correct model category
      std::smatch sm;
      std::regex_search(fullParamName, sm, modelPattern);
      std::string modelSuffix = sm[sm.size()-1];
      uint modelId = static_cast<uint>(stoi(modelSuffix));
      if (mapSharedParams){
        int type = getEnumOfParamName(pattern);
        auto it = mapSharedParams->begin();
        while (it != mapSharedParams->end()){
            auto paramNum = it->first;
            auto sharedParams = (*mapSharedParams)[paramNum];
            std::pair<uint, int> modelAndType;
            modelAndType.first = modelId;
            modelAndType.second = type;
            auto found = std::find(sharedParams.begin(), sharedParams.end(), modelAndType);
            if (found != sharedParams.end()){
                for (size_t j = 0; j < (*mapSharedParams)[paramNum].size(); j ++){
                    uint model_j = (*mapSharedParams)[paramNum][j].first;
                    matchingParamsPerModel[model_j].push_back(fullParamName);
                }
                break;
            }
            it ++;
        }
      }
    }
  }

  return matchingParamsPerModel;

}

// /*******************************************************************************/
std::vector <string> ChromosomeNumberOptimizer::getNonFixedParams(SingleProcessPhyloLikelihood* tl, ParameterList &allParams, std::map<uint, vector<int>>* fixedParams) const{
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, tl);
    uint numOfModels = static_cast<uint>(tl->getSubstitutionProcess().getNumberOfModels());
    vector<string> nonFixed;
    for (int i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
        auto it = typeWithParamNames.find(i);
        if (it == typeWithParamNames.end()){
            continue;
        }
        int type = it->first;
        auto modelAndParameterNames = typeWithParamNames[type];
        for(uint j = 1; j <= numOfModels; j ++){
            if (!(std::count((*fixedParams)[j].begin(), (*fixedParams)[j].end(), type))){
                vector<string> parameterNames = modelAndParameterNames[j];
                for (size_t k = 0; k < parameterNames.size(); k++){
                    nonFixed.push_back(parameterNames[k]);
                }
            }

        }
        
    }
    return nonFixed;
}

// /***************************************************************************************/
void ChromosomeNumberOptimizer::optimizeBaseNum(SingleProcessPhyloLikelihood* tl, size_t index, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, double lowerBound, 
                                                double upperBound, const string &paramName, ParameterList& params, uint model, std::map<uint, uint> &baseNumberUpperBounds){

    Function* func = tl;
    ParameterList substitutionParams = tl->getSubstitutionModelParameters();
    vector <std::string> names;
    for (size_t k = 0; k < substitutionParams.size(); k ++){
        names.push_back(substitutionParams[k].getName());
    }
    ParameterList updatedSubstitutionParams = params.createSubList(names);
    size_t best_i = (size_t)(params.getParameter(paramName).getValue());
    double f_value = *currentLikelihood;
    
    for (size_t i = 0; i < baseNumCandidates.size(); i++){
        unsigned int baseNum = baseNumCandidates[i];
        if (baseNum > baseNumberUpperBounds[model]){
            break;
        }
        params.getParameter(paramName).setValue((double)baseNum);
        double f_i = func->f(params);
        if (f_i < f_value){
            best_i = baseNum;
            f_value = f_i;
        }
    }
    params.getParameter(paramName).setValue((double)best_i);
    updatedSubstitutionParams.getParameter(paramName).setValue((double)best_i);
    func->f(updatedSubstitutionParams);
    //param.setValue((double)best_i);
    *currentLikelihood = f_value;

}
/******************************************
Functions for heterogeneous ChromEvol model
*******************************************/
double ChromosomeNumberOptimizer::calculateModelSelectionCriterion(SingleProcessPhyloLikelihood* lik, size_t numOfFixedParams) const{
    double modelSelectionCriterionValue;
    if (ChromEvolOptions::modelSelectionCriterion_ == "AICc"){
        modelSelectionCriterionValue = calculateAICc(lik, numOfFixedParams);

    }else if (ChromEvolOptions::modelSelectionCriterion_ == "AIC"){
        modelSelectionCriterionValue = calculateAIC(lik, numOfFixedParams);

    }else{
        throw Exception("ChromosomeNumberOptimizer::calculateModelSelectionCriterion(): No such criterion exists!");
    }
    return modelSelectionCriterionValue;
}

double ChromosomeNumberOptimizer::calculateAICc(SingleProcessPhyloLikelihood* lik, size_t numOfFixedParams) const{
    // the number of shifts
    size_t numOfModels = lik->getSubstitutionProcess().getNumberOfModels();
    // the number of substitution params takes into account also the backward phase
    auto numOfSubstitutionParams = lik->getSubstitutionModelParameters().size()-numOfFixedParams;

    // N (sample size)
    auto sampleSize = tree_->getAllLeavesNames().size();
    // p (number of overall parameters)
    // k must maintain the constraint of the denominator -> throw detailed exception..
    double numOfParams = static_cast<double>(numOfModels) - 1 + static_cast<double>(numOfSubstitutionParams);
    if ((backwardPhaseStarted_) && (numOfShiftsForward_ > numOfModels)){
        numOfParams += static_cast<double>(numOfShiftsForward_-numOfModels);
    }
    // sample size correction term
    auto denominator = (static_cast<double>(sampleSize)-numOfParams-1);
    if (denominator <= 0){
        string errorMessage = "ChromosomeNumberOptimizer::calculateAICc(): Too many parameters when using " + std::to_string(numOfModels) + " different models! Try to set _maxNumOfModels to " + std::to_string(numOfModels-1) + "\n"; 
        throw Exception(errorMessage);
    }
    double sampleSizeCorrection = (2*numOfParams*(numOfParams + 1))/denominator;
    //Calculating AICc-> I have some problems with AIC. In some cases the denominator becomes zero.
    // so meanwhile I will use AIC...
    double AIC = 2*(lik->getValue()) + (2*numOfParams);
    double AICc = AIC + sampleSizeCorrection;
    //return AICc;
    return AICc;

}
double ChromosomeNumberOptimizer::calculateAIC(SingleProcessPhyloLikelihood* lik, size_t numOfFixedParams) const{
    // the number of shifts
    size_t numOfModels = lik->getSubstitutionProcess().getNumberOfModels();
    // the number of substitution params takes into account also the backward phase
    auto numOfSubstitutionParams = lik->getSubstitutionModelParameters().size()-numOfFixedParams;

    // N (sample size)
    // p (number of overall parameters)
    // k must maintain the constraint of the denominator -> throw detailed exception..
    double numOfParams = static_cast<double>(numOfModels) - 1 + static_cast<double>(numOfSubstitutionParams);
    if ((backwardPhaseStarted_) && (numOfShiftsForward_ > numOfModels)){
        numOfParams += static_cast<double>(numOfShiftsForward_-numOfModels);
    }

    //Calculating AICc-> I have some problems with AIC. In some cases the denominator becomes zero.
    // so meanwhile I will use AIC...
    double AIC = 2*(lik->getValue()) + (2*numOfParams);
    //return AICc;
    return AIC;

}
/***********************************************/
void ChromosomeNumberOptimizer::setInitialModelRepresentitives(std::map<uint, vector<uint>> &initialPartition){
    auto numOfModels = static_cast<uint>(initialPartition.size());
    for (uint i = 0; i < numOfModels; i++){
        prevModelsPartitions_[numOfModels].push_back(initialPartition[i+1][0]);
    }
}
/***********************************************/
void ChromosomeNumberOptimizer::setParamsNameInForMultiProcess(std::map<uint, std::map<int, vector<string>>> &mapOfParamsNamesPerModelType, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams){
    auto it = modelParams.begin();
    while(it != modelParams.end()){
        uint model = it->first;
        for (int i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
            vector<string> paramsFullNames;
            auto strParamName = getStringParamName(i);
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
/***********************************************/
void ChromosomeNumberOptimizer::setNewModelAttributes(SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, std::map<uint, pair<int, std::map<int, std::vector<double>>>>* modelParams, std::map<uint, std::vector<uint>>* mapModelNodesIds, std::map<uint, uint>* baseNumberBounds){
    uint modelNumInCurrentLik = static_cast<uint>(currentLik->getSubstitutionProcess().getModelNumberForNode(nodeToSplit));
    uint numOfModels = static_cast<uint>(currentLik->getSubstitutionProcess().getNumberOfModels());
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, 0, currentLik, sharedParams);   
    *modelParams = getMapOfParamsForComplexModel(currentLik, typeWithParamNames, numOfModels);

    // add the new regime to the map (i.e., copy the parameters of the regime to which the node was assigned to previously)
    (*modelParams)[numOfModels+1] = (*modelParams)[modelNumInCurrentLik];
    //getting nodeIds per each model (for the current model) *not the new one
    getMutableMapOfModelAndNodeIds(*mapModelNodesIds, currentLik);

    // getting an updated map of nodes per model (for the new model)
    PhyloTree tree = *tree_;
    vector<std::shared_ptr<PhyloNode>> newSubtree = tree.getSubtreeNodes(tree.getNode(nodeToSplit));
    vector<uint> newSubtreeIds = tree.getNodeIndexes(newSubtree);
    //delete tree;
    (*mapModelNodesIds)[numOfModels+1] = std::vector<uint>();
    for (size_t i = 0; i < newSubtreeIds.size(); i++){
        uint nodeId = newSubtreeIds[i];
        uint prevModel = static_cast<uint>(currentLik->getSubstitutionProcess().getModelNumberForNode(nodeId));
        if (prevModel == modelNumInCurrentLik){
            (*mapModelNodesIds)[prevModel].erase(std::remove((*mapModelNodesIds)[prevModel].begin(), (*mapModelNodesIds)[prevModel].end(), nodeId), (*mapModelNodesIds)[prevModel].end());
            (*mapModelNodesIds)[numOfModels+1].push_back(nodeId);

        }
    }
    for (size_t i = 1; i <= numOfModels; i ++){
        auto branchProcess = currentLik->getSubstitutionProcess().getModel(i);
        baseNumberUpperBound_[static_cast<uint>(i)] = std::dynamic_pointer_cast<const ChromosomeSubstitutionModel>(branchProcess)->getMaxChrRange();
        (*baseNumberBounds)[static_cast<uint>(i)] = baseNumberUpperBound_[static_cast<uint>(i)];
    }
    baseNumberUpperBound_[numOfModels+1] = baseNumberUpperBound_[modelNumInCurrentLik];
    (*baseNumberBounds)[numOfModels+1] = baseNumberUpperBound_[modelNumInCurrentLik];

}
/***********************************************/

void ChromosomeNumberOptimizer::optimizeFirstRound(std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, double parsimonyBound, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, vector<uint> numOfPointsNextRounds, vector<uint> numOfIterationsNextRounds, vector<SingleProcessPhyloLikelihood*> &vectorOfLikelihoods, string* text, std::map<uint, uint>* baseNumberBounds, std::map<uint, uint>* mapOfModelsBackward, std::map<uint, pair<int, std::map<int, std::vector<double>>>>* prevModelParamsBackward, std::map<uint, vector<uint>>* modelsBackwards, omp_lock_t* mutex){
    size_t index = min((int)numOfPointsNextRounds.size()-1, 1);
    vectorOfLikelihoods.reserve(numOfPointsNextRounds[index]);
    if (mutex){
        omp_set_lock(mutex);
    } 
    vector <unsigned int> baseNumCandidates;
    getBaseNumCandidates(baseNumCandidates, *baseNumberBounds);
    if (mutex){
        omp_unset_lock(mutex);

    }  
    string log = "";

    // If base number is one of the parameters
    log += "##################################\n";
    log += "*********  cycle 0  **************\n";
    printLog(text, log);  
    for (size_t n = 0; n < numOfPointsNextRounds[0]; n++){
        log = "Starting cycle with Point #";
        log += std::to_string(n) +"....\n";
        printLog(text, log);
        if (mutex){
            omp_set_lock(mutex);
        }    //std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, uint numOfPoints, std::map<uint, vector<int>> &fixedParams, double parsimonyBound, uint iteration
        //std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, double parsimonyBound, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, uint iteration
        SingleProcessPhyloLikelihood* lik;
        if (!(mapOfModelsBackward)){
            lik = getSingleNewLikObject(updatedSharedParams, fixedParams, parsimonyBound, mapModelNodesIds, modelParams, numOfModels, (uint)n, baseNumberBounds);
        }else{
            lik = getBackwardLikObject(updatedSharedParams, fixedParams, parsimonyBound, mapModelNodesIds, modelParams, numOfModels, *prevModelParamsBackward, *mapOfModelsBackward, *modelsBackwards, (uint)n, baseNumberBounds);

        }
        if (mutex){
            omp_unset_lock(mutex);

        }      
        fillVectorOfLikelihoods(lik, numOfIterationsNextRounds[0],  n, numOfPointsNextRounds[index], baseNumCandidates, updatedSharedParams, fixedParams, vectorOfLikelihoods, text, *baseNumberBounds, mutex);
   
    }
    if (mutex){
        omp_set_lock(mutex);
    } 
    sort(vectorOfLikelihoods.begin(), vectorOfLikelihoods.end(), compareLikValues);
    if (mutex){
        omp_unset_lock(mutex);

    }
    printLikelihoodVectorValues(vectorOfLikelihoods, text, 0);

}
/***********************************************/
size_t ChromosomeNumberOptimizer::getMaxNumOfMergingModels(std::map<uint, vector<uint>> &modelsToMerge){
    size_t maxNumOfModels = 0;
    auto it = modelsToMerge.begin();
    while (it != modelsToMerge.end()){
        auto sizeOfMergedCluster = modelsToMerge[it->first].size();
        if (sizeOfMergedCluster == 0){
            it ++;
            continue;
        }
        sizeOfMergedCluster += 1;
        if (sizeOfMergedCluster > maxNumOfModels){
            maxNumOfModels = sizeOfMergedCluster;
        }
 
        it ++;
    }
    return maxNumOfModels;
}
/***********************************************/
SingleProcessPhyloLikelihood* ChromosomeNumberOptimizer::getBackwardLikObject(std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, double parsimonyBound, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &prevModelParams, std::map<uint, uint> &mapOfModelsToMerge, std::map<uint, vector<uint>> &modelsToMerge, uint iteration, std::map<uint, uint>* baseNumberBounds){
    SingleProcessPhyloLikelihood* newLik;
    size_t maxNumOfMergingModels = getMaxNumOfMergingModels(modelsToMerge);
    if (iteration <= maxNumOfMergingModels){
        auto it = modelsToMerge.begin();
        while(it != modelsToMerge.end()){
            if (modelsToMerge[it->first].size() == 0){
                it ++;
                continue;
            }
            if (iteration < modelsToMerge[it->first].size() + 1){
                if (iteration == modelsToMerge[it->first].size()){
                    modelParams[mapOfModelsToMerge[it->first]] = prevModelParams[it->first];  
                }else{
                    modelParams[mapOfModelsToMerge[it->first]] = prevModelParams[modelsToMerge[it->first][iteration]];

                }
                
            }else{
                // get mean for each parameter
                modelParams[mapOfModelsToMerge[it->first]] = getMeanParameters(prevModelParams, modelsToMerge[it->first], it->first);
            }

            it ++;
        }

        newLik = setHeterogeneousModel(tree_, vsc_, alphabet_, *baseNumberBounds, mapModelNodesIds, modelParams, numOfModels, updatedSharedParams);         

    }else{
        if (iteration == 1){
            newLik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, *baseNumberBounds, mapModelNodesIds, modelParams, numOfModels, parsimonyBound * (double)iteration, fixedParams, updatedSharedParams);

        }else{
            newLik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, *baseNumberBounds, mapModelNodesIds, modelParams, numOfModels, parsimonyBound * (1+ (0.1*(double)iteration)), fixedParams, updatedSharedParams);
        }
        
    }
    ifNanTryToResampleLikObject(&newLik, tree_, vsc_, alphabet_, *baseNumberBounds, mapModelNodesIds, modelParams, numOfModels, parsimonyBound, iteration, fixedParams, updatedSharedParams);
    return newLik;

}
/***********************************************/
SingleProcessPhyloLikelihood* ChromosomeNumberOptimizer::getSingleNewLikObject(std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, double parsimonyBound, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, uint iteration, std::map<uint, uint>* baseNumberBounds){

    SingleProcessPhyloLikelihood* newLik;
    if (iteration == 0){
        newLik = setHeterogeneousModel(tree_, vsc_, alphabet_, *baseNumberBounds, mapModelNodesIds, modelParams, numOfModels+1, updatedSharedParams);         

    }else{
        if (iteration == 1){
            newLik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, *baseNumberBounds, mapModelNodesIds, modelParams, numOfModels+1, parsimonyBound * (double)iteration, fixedParams, updatedSharedParams);

        }else{
            newLik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, *baseNumberBounds, mapModelNodesIds, modelParams, numOfModels+1, parsimonyBound * (1+(0.1*(double)iteration)), fixedParams, updatedSharedParams);
        }
        
    }
    ifNanTryToResampleLikObject(&newLik, tree_, vsc_, alphabet_, *baseNumberBounds, mapModelNodesIds, modelParams, numOfModels+1, parsimonyBound, iteration, fixedParams, updatedSharedParams);
    return newLik;
    

}
/***********************************************/
// this function will replace the getNewLikObject() once it is tested for bugs
void ChromosomeNumberOptimizer::getNewLikObjectForParallelRuns(std::vector<SingleProcessPhyloLikelihood*> &perCandidateLikVec, SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, uint numOfPoints, std::map<uint, vector<int>> &fixedParams, double parsimonyBound){
    
    uint modelNumInCurrentLik = static_cast<uint>(currentLik->getSubstitutionProcess().getModelNumberForNode(nodeToSplit));
    uint numOfModels = static_cast<uint>(currentLik->getSubstitutionProcess().getNumberOfModels());
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, 0, currentLik, sharedParams);
    
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> modelParams = getMapOfParamsForComplexModel(currentLik, typeWithParamNames, numOfModels);

    // add the new regime to the map (i.e., copy the parameters of the regime to which the node was assigned to previously)
    modelParams[numOfModels+1] = modelParams[modelNumInCurrentLik];
    //getting nodeIds per each model (for the current model) *not the new one
    std::map<uint, std::vector<uint>> mapModelNodesIds;
    getMutableMapOfModelAndNodeIds(mapModelNodesIds, currentLik);

    // getting an updated map of nodes per model (for the new model)
    PhyloTree tree = *tree_;
    vector<std::shared_ptr<PhyloNode>> newSubtree = tree.getSubtreeNodes(tree.getNode(nodeToSplit));
    vector<uint> newSubtreeIds = tree.getNodeIndexes(newSubtree);
    //delete tree;
    mapModelNodesIds[numOfModels+1] = std::vector<uint>();
    for (size_t i = 0; i < newSubtreeIds.size(); i++){
        uint nodeId = newSubtreeIds[i];
        uint prevModel = static_cast<uint>(currentLik->getSubstitutionProcess().getModelNumberForNode(nodeId));
        if (prevModel == modelNumInCurrentLik){
            mapModelNodesIds[prevModel].erase(std::remove(mapModelNodesIds[prevModel].begin(), mapModelNodesIds[prevModel].end(), nodeId), mapModelNodesIds[prevModel].end());
            mapModelNodesIds[numOfModels+1].push_back(nodeId);

        }
    }
    for (size_t i = 1; i <= numOfModels; i ++){
        auto branchProcess = currentLik->getSubstitutionProcess().getModel(i);
        baseNumberUpperBound_[static_cast<uint>(i)] = std::dynamic_pointer_cast<const ChromosomeSubstitutionModel>(branchProcess)->getMaxChrRange();
    }
    baseNumberUpperBound_[numOfModels+1] = baseNumberUpperBound_[modelNumInCurrentLik];
    // setting the heterogeneous model
    SingleProcessPhyloLikelihood* newLik;
    for (size_t i = 0; i < numOfPoints; i++){
        if (i == 0){
            newLik = setHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels+1, updatedSharedParams);

        }else{
            if (i == 1){
                newLik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels+1, parsimonyBound * (double)i, fixedParams, updatedSharedParams);

            }else{
                newLik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels+1, parsimonyBound * (1+ ((double)i*0.1)), fixedParams, updatedSharedParams);
            }
            
        }
        ifNanTryToResampleLikObject(&newLik, tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels+1, parsimonyBound, numOfPoints, fixedParams, updatedSharedParams);
        perCandidateLikVec.push_back(newLik);
    }
    
}
/***********************************************/

void ChromosomeNumberOptimizer::getNewLikObject(SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, uint numOfPoints, std::map<uint, vector<int>> &fixedParams, double parsimonyBound){
    
    uint modelNumInCurrentLik = static_cast<uint>(currentLik->getSubstitutionProcess().getModelNumberForNode(nodeToSplit));
    uint numOfModels = static_cast<uint>(currentLik->getSubstitutionProcess().getNumberOfModels());
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, 0, currentLik, sharedParams);
    
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> modelParams = getMapOfParamsForComplexModel(currentLik, typeWithParamNames, numOfModels);

    // add the new regime to the map (i.e., copy the parameters of the regime to which the node was assigned to previously)
    modelParams[numOfModels+1] = modelParams[modelNumInCurrentLik];
    //getting nodeIds per each model (for the current model) *not the new one
    std::map<uint, std::vector<uint>> mapModelNodesIds;
    getMutableMapOfModelAndNodeIds(mapModelNodesIds, currentLik);

    // getting an updated map of nodes per model (for the new model)
    PhyloTree tree = *tree_;
    vector<std::shared_ptr<PhyloNode>> newSubtree = tree.getSubtreeNodes(tree.getNode(nodeToSplit));
    vector<uint> newSubtreeIds = tree.getNodeIndexes(newSubtree);
    mapModelNodesIds[numOfModels+1] = std::vector<uint>();
    for (size_t i = 0; i < newSubtreeIds.size(); i++){
        uint nodeId = newSubtreeIds[i];
        uint prevModel = static_cast<uint>(currentLik->getSubstitutionProcess().getModelNumberForNode(nodeId));
        if (prevModel == modelNumInCurrentLik){
            mapModelNodesIds[prevModel].erase(std::remove(mapModelNodesIds[prevModel].begin(), mapModelNodesIds[prevModel].end(), nodeId), mapModelNodesIds[prevModel].end());
            mapModelNodesIds[numOfModels+1].push_back(nodeId);

        }
    }
    for (size_t i = 1; i <= numOfModels; i ++){
        auto branchProcess = currentLik->getSubstitutionProcess().getModel(i);
        baseNumberUpperBound_[static_cast<uint>(i)] = std::dynamic_pointer_cast<const ChromosomeSubstitutionModel>(branchProcess)->getMaxChrRange();
    }
    baseNumberUpperBound_[numOfModels+1] = baseNumberUpperBound_[modelNumInCurrentLik];
    // setting the heterogeneous model
    SingleProcessPhyloLikelihood* newLik;
    for (size_t i = 0; i < numOfPoints; i++){
        if (i == 0){
            newLik = setHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels+1, updatedSharedParams);

        }else{
            if (i == 1){
                newLik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels+1, parsimonyBound * (double)i, fixedParams, updatedSharedParams);

            }else{
                newLik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels+1, parsimonyBound * (1+(0.1*(double)i)), fixedParams, updatedSharedParams);
            }
            
        }
        ifNanTryToResampleLikObject(&newLik, tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels+1, parsimonyBound, numOfPoints, fixedParams, updatedSharedParams);
        vectorOfLikelohoods_.push_back(newLik);
    }
    
}
/****************************************************************/
std::map<uint, pair<int, std::map<int, std::vector<double>>>> ChromosomeNumberOptimizer::getMapOfParamsForComplexModel(SingleProcessPhyloLikelihood* lik, std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames, uint numOfModels) {
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
/***********************************************************************************************/
//ifNanTryToResampleLikObject(lik, tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels, parsimonyBound, numOfPoints, fixedParams, sharedParams);
void ChromosomeNumberOptimizer::ifNanTryToResampleLikObject(SingleProcessPhyloLikelihood** lik, const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet,
    std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, 
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, 
    uint numOfModels, double parsimonyBound, int numOfPoints,
    std::map<uint, vector<int>> &fixedParams,
    std::map<int, std::vector<std::pair<uint, int>>>* sharedParams)
{
    int numOfTrials = 0;
    while((std::isnan((*lik)->getValue())) && (numOfTrials < ChromEvolOptions::maxNumOfTrials_)){
        std::cerr << "WARNING!!! The value of the likelihood object is nan"<< std::endl;
        auto likToDel = *lik; 
        deleteLikObject(likToDel);
        int factor = 0;
        if (numOfPoints == 0){
            factor ++;
        }else{
            factor = numOfPoints;
        }
        auto multiplier = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(factor);
        if (multiplier == 0){
            multiplier ++;
        }
        double parsimonlyFactor;
        if (multiplier == 1){
            parsimonlyFactor = (double)multiplier;
        }else{
            parsimonlyFactor = (1 + (0.1 * (double)multiplier));
        }
        *lik = setRandomHeterogeneousModel(tree, vsc, alphabet, baseNumberUpperBound, mapModelNodesIds, modelParams, numOfModels, parsimonyBound * parsimonlyFactor, fixedParams, sharedParams);
        numOfTrials ++;
    }
    if (std::isnan((*lik)->getValue())){
        throw Exception("ERROR!!! The likelihood point is nan!!!!\n");
    }
}

/***********************************************************************************************/
void ChromosomeNumberOptimizer::getMutableMapOfModelAndNodeIds(std::map<uint, vector<uint>> &mapModelNodesIds, SingleProcessPhyloLikelihood* lik, uint rootId){
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
/**********************************************************************************************/
SingleProcessPhyloLikelihood* ChromosomeNumberOptimizer::setRandomHeterogeneousModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet,
    std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, 
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, 
    uint numOfModels, double parsimonyBound, 
    std::map<uint, vector<int>> &fixedParams,
    std::map<int, std::vector<std::pair<uint, int>>>* sharedParams)
{


    std::map<uint, std::map<int, vector<string>>> mapOfParamsNamesPerModelType;
    setParamsNameInForMultiProcess(mapOfParamsNamesPerModelType, modelParams);
    std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
    string fixedRootFreqPath = ChromEvolOptions::fixedFrequenciesFilePath_;
    bool weightedRootFreqs;
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::shared_ptr<ChromosomeSubstitutionModel>(ChromosomeSubstitutionModel::initRandomModel(alphabet, modelParams[1].first, modelParams[1].second, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_, fixedParams[1], parsimonyBound));
    if (fixedRootFreqPath == "none"){
        weightedRootFreqs = true;
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree);

    }else{
        weightedRootFreqs = false;
        vector <double> rootFreqs = setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
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
    }
    aliasParametersInSubstitutionProcess(mapOfParamsNamesPerModelType, sharedParams, subProSim);
    SubstitutionProcess* nsubPro= subProSim->clone();
    Context* context = new Context();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;

}
/**********************************************************************************************/
void ChromosomeNumberOptimizer::aliasParametersInSubstitutionProcess(std::map<uint, std::map<int, vector<string>>> &mapOfParamsNamesPerModelType, std::map<int, vector<std::pair<uint, int>>>* updatedSharedParams, std::shared_ptr<NonHomogeneousSubstitutionProcess> process){
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
/**********************************************************************************************/
SingleProcessPhyloLikelihood* ChromosomeNumberOptimizer::setHeterogeneousModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, std::map<int, vector<std::pair<uint, int>>>* updatedSharedParams){
    std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
    string fixedRootFreqPath = ChromEvolOptions::fixedFrequenciesFilePath_;
    bool weightedRootFreqs;
    std::map<uint, std::map<int, vector<string>>> mapOfParamsNamesPerModelType;
    setParamsNameInForMultiProcess(mapOfParamsNamesPerModelType, modelParams);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, modelParams[1].second, modelParams[1].first, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
    if (fixedRootFreqPath == "none"){
        weightedRootFreqs = true;
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree);

    }else{
        weightedRootFreqs = false;
        vector <double> rootFreqs = setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
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
    }
    aliasParametersInSubstitutionProcess(mapOfParamsNamesPerModelType, updatedSharedParams, subProSim);


    SubstitutionProcess* nsubPro= subProSim->clone();
    Context* context = new Context();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;

}
/**********************************************************************************************/
void ChromosomeNumberOptimizer::updateSharedParameters(std::map<int, vector<std::pair<uint, int>>> &sharedParams, uint prevShift, uint numOfShifts) const{
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
                //int paramNumToPut;
                if (isSharedBetweenModels){
                // if for example gain2 = gain1, if model 3 derives from model2, gain3 = gain1.
                    sharedParams[paramNum].push_back(modelAndParam);
                    //paramNumToPut = paramNum;
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
/**********************************************************************************************/

// this function will replace the optimize function once it will be tested on enough data
void ChromosomeNumberOptimizer::optimizeInParallel(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds){
    int maxNumOfModels;
    SingleProcessPhyloLikelihood* lik;
    SingleProcessPhyloLikelihood* minAICcLik;
    
    fixedParams_ = fixedParams;
    sharedParams_ = ChromEvolOptions::sharedParameters_;
    optimizeBaseNumber_ = false;
    // if we should optimize it for at least one model, we will set it to true. Otherwise it is fixed for all the models.
    for (uint i = 1; i <= static_cast<uint>(ChromEvolOptions::numOfModels_); i++){
        if (!(std::count(fixedParams_[i].begin(), fixedParams_[i].end(), ChromosomeSubstitutionModel::BASENUM))){
            optimizeBaseNumber_ = true;
        }
    }
    ((ChromEvolOptions::maxNumOfModels_ == 1) && (ChromEvolOptions::heterogeneousModel_)) ? (maxNumOfModels = (static_cast<int>((tree_->getAllLeavesNames()).size())-1)) : (maxNumOfModels = ChromEvolOptions::maxNumOfModels_);
    vector<uint> candidateShiftNodesIds;
    if (ChromEvolOptions::heterogeneousModel_){
        getValidCandidatesForShift(candidateShiftNodesIds, ChromEvolOptions::minCladeSize_, ChromEvolOptions::numOfModels_);

    }
    
    vector <unsigned int> baseNumCandidates;
    getBaseNumCandidates(baseNumCandidates, baseNumberUpperBound_);
    uint numOfShifts = ChromEvolOptions::numOfModels_;

    sharedParams_ = ChromEvolOptions::sharedParameters_;
    vectorOfLikelohoods_.reserve(numOfPoints);
    //vectorOfContexts_.reserve(numberOfModels);

    if (seed != 0){
        RandomTools::setSeed(static_cast<long>(seed));
    }
    bool deltaAICcImproved = true;
    bool firstIteration = true;

    while((deltaAICcImproved) && (numOfShifts <= (size_t)maxNumOfModels)){
        if ((candidateShiftNodesIds.size() == 0) && (maxNumOfModels > 1)){
            break;
        }
        if (firstIteration){
            initLikelihoods(modelParams, parsimonyBound, rateChange, numOfPoints, fixedRootFreqPath, fixedParams_, mapModelNodesIds, numOfShifts, &sharedParams_);
            optimizeMultiProcessModel(&sharedParams_, &fixedParams_, numOfPoints_, numOfIterations_, baseNumberUpperBound_, 0,0);
            // leave only the best one
            clearVectorOfLikelihoods(1);
            firstIteration = false;
            if ((prevModelsAICcLikValues_.empty()) && (ChromEvolOptions::heterogeneousModel_)){
                size_t numOfFixedParams = getNumberOfFixedParams(vectorOfLikelohoods_[0], fixedParams_); 
                double AICc = calculateModelSelectionCriterion(vectorOfLikelohoods_[0], numOfFixedParams);
                std::pair<double, double> AICcAndLik(AICc, vectorOfLikelohoods_[0]->getValue());
                prevModelsAICcLikValues_[ChromEvolOptions::numOfModels_] = AICcAndLik;
                getParameterNamesAndValues(vectorOfLikelohoods_[0], ChromEvolOptions::numOfModels_);
                prevModelsRootFrequencies_[ChromEvolOptions::numOfModels_] = getRootFrequencies(vectorOfLikelohoods_[0]);
            }
            numOfShifts ++;
            continue;

        }
        
        lik = vectorOfLikelohoods_.back();
        minAICcLik = lik;
        auto prevLik = minAICcLik;

        //vectorOfLikelohoods_.pop_back();     
        size_t numOfFixedParams = getNumberOfFixedParams(minAICcLik, fixedParams_); 
        double initialAICc = calculateModelSelectionCriterion(minAICcLik, numOfFixedParams);
        uint minDetaAICcNode;
        std::cout << "*** *** *** Starting considering " << numOfShifts << " shifts *** *** ***" << std::endl;
        //omp_set_num_threads(4);
        SingleProcessPhyloLikelihood* bestLikAmongCandidates = 0;
        double bestAICcScore = initialAICc;
        omp_lock_t mutex;
        omp_init_lock(&mutex);
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < candidateShiftNodesIds.size(); i++){
            runNewBranchModel(mutex, lik, candidateShiftNodesIds, i, numOfShifts, parsimonyBound, numOfPoints, &bestLikAmongCandidates, &bestAICcScore, &minDetaAICcNode);
            std::cout << "Number of threads in iteration  " << i << " : " << omp_get_num_threads() << std::endl;
        }
        omp_destroy_lock(&mutex);


        // get the best candidate
        // for each candidate calculate the AICc, and delete those which have a worse AICc than the best so far
        uint prevShiftOfBest = static_cast<uint>(lik->getSubstitutionProcess().getModelNumberForNode(minDetaAICcNode));
        std::map<uint, vector<int>> fixedParameters = fixedParams_;
        if (initialAICc - bestAICcScore > ChromEvolOptions::deltaAICcThreshold_){
            //minAICc = bestAICcScore;
            std::cout << "***" << std::endl;
            vectorOfLikelohoods_.pop_back();
            deleteLikObject(prevLik);
            minAICcLik = bestLikAmongCandidates;
            deltaAICcImproved = true;
            std::map<int, std::vector<std::pair<uint, int>>> sharedParams = sharedParams_;
            updateSharedParameters(sharedParams, prevShiftOfBest, numOfShifts);
            sharedParams_ = sharedParams;
            fixedParameters[numOfShifts] = fixedParameters[prevShiftOfBest];
            prevModelsPartitions_[numOfShifts] = prevModelsPartitions_[numOfShifts-1];
            prevModelsPartitions_[numOfShifts].push_back(minDetaAICcNode);
            std::pair<double, double> AICcAndLik(bestAICcScore, minAICcLik->getValue());
            prevModelsAICcLikValues_[numOfShifts] = AICcAndLik;
            prevModelsRootFrequencies_[numOfShifts] = getRootFrequencies(minAICcLik);
            getParameterNamesAndValues(minAICcLik, numOfShifts);
            fixedParams_ = fixedParameters;
            
            vectorOfLikelohoods_.push_back(minAICcLik);
            candidateShiftNodesIds.clear();
            getValidCandidatesForShift(candidateShiftNodesIds, ChromEvolOptions::minCladeSize_, numOfShifts);
            numOfShifts ++;
            //candidateShiftNodesIds.erase(std::remove(candidateShiftNodesIds.begin(), candidateShiftNodesIds.end(), minDetaAICcNode), candidateShiftNodesIds.end());
        }else{
            deleteLikObject(bestLikAmongCandidates);
            deltaAICcImproved = false;

        }
        
        
    }
     std::cout << "*** Final best model: " << vectorOfLikelohoods_[0]->getValue() << std::endl;

    //initLikelihoods(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds, uint numOfModels, std::map<int, vector<uint>>* sharedParams, std::map<uint, vector<int>> &fixedParameters)

}
/**********************************************************************************************/
void ChromosomeNumberOptimizer::getParameterNamesAndValues(SingleProcessPhyloLikelihood* lik, uint numOfModels){
     auto params = lik->getSubstitutionModelParameters();
     for (size_t i = 0; i < params.size(); i++){
         auto paramName = params[i].getName();
         double value = lik->getLikelihoodCalculation()->getParameter(paramName).getValue();
         std::pair<string, double> paramAndValue(paramName, value);
         prevModelParams_[numOfModels].push_back(paramAndValue);

     }
}
/**********************************************************************************************/
void ChromosomeNumberOptimizer::runNewBranchModel(omp_lock_t &mutex, SingleProcessPhyloLikelihood* lik, vector<uint> &candidateShiftNodesIds, size_t i, uint numOfShifts, double parsimonyBound, uint numOfPoints, SingleProcessPhyloLikelihood** bestCandidateLik, double* bestAICc, uint* minAICcNode){
    vector<SingleProcessPhyloLikelihood*> perCandidateLikVec; // this vector should hold the likelihoods related to the cucles of the currently examined candidate
    std::map<int, std::vector<std::pair<uint, int>>> sharedParams = sharedParams_;
    std::map<uint, vector<int>> fixedParams = fixedParams_;
    uint prevShift = static_cast<uint>(lik->getSubstitutionProcess().getModelNumberForNode(candidateShiftNodesIds[i]));
    updateSharedParameters(sharedParams, prevShift, numOfShifts);
    fixedParams[numOfShifts] = fixedParams[prevShift];
    // the following section is a critical section, because the tree_ object adds to the
    // observers_ data member the parametrizable tree. Without this mutex, there are segmentation faults.
    
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> modelParams;
    std::map<uint, std::vector<uint>> mapModelNodesIds;
    omp_set_lock(&mutex);
    //std::cout << "Within critical section before setNewModelAttributes()" << std::endl;
    std::map<uint, uint> baseNumberBounds;
    setNewModelAttributes(lik, candidateShiftNodesIds[i], &sharedParams_, &sharedParams, fixedParams, &modelParams, &mapModelNodesIds, &baseNumberBounds);
    //std::cout << "Out of critical section before setNewModelAttributes()" << std::endl;
    omp_unset_lock(&mutex);
    uint numOfModels =  static_cast<uint>(lik->getSubstitutionProcess().getNumberOfModels());
    string textToPrint = "";
    optimizeFirstRound(&sharedParams, fixedParams, parsimonyBound, mapModelNodesIds, modelParams, numOfModels, numOfPointsNextRounds_, numOfIterationsNextRounds_, perCandidateLikVec, &textToPrint, &baseNumberBounds, 0, 0, 0, &mutex);

    optimizeMultiProcessModel(&sharedParams, &fixedParams, numOfPointsNextRounds_, numOfIterationsNextRounds_, baseNumberBounds, &perCandidateLikVec, &textToPrint, &mutex);
    //std::cout << "After optimizeMultiProcessModel: " << i << std::endl;
    size_t numOfFixedParams = getNumberOfFixedParams(perCandidateLikVec[0], fixedParams); 
    double AICc = calculateModelSelectionCriterion(perCandidateLikVec[0], numOfFixedParams);
    omp_set_lock(&mutex);
    std::cout << "\tshift is at node: "<< candidateShiftNodesIds[i] << std::endl;
    std::cout << textToPrint << std::endl;
    std::cout << "Final log likelihood is: " << perCandidateLikVec[0]->getValue() << std::endl;
    std::cout << "Final "<< ChromEvolOptions::modelSelectionCriterion_<< " is: " << AICc << std::endl;
    std::cout << "*** *** *** *** "<< candidateShiftNodesIds[i] << std::endl;
    printLikParameters(perCandidateLikVec[0], 1, &textToPrint);
    if ((*bestCandidateLik == 0) || (AICc < *bestAICc)){                      
        *bestAICc =  AICc;
        *minAICcNode = candidateShiftNodesIds[i];
        if (*bestCandidateLik != 0){
            auto likObjectToDel = *bestCandidateLik;
            *bestCandidateLik = perCandidateLikVec[0];
            deleteLikObject(likObjectToDel);
        }else{
            *bestCandidateLik = perCandidateLikVec[0];
        }

    }else{
        SingleProcessPhyloLikelihood* lik_to_del = perCandidateLikVec.back(); 
        perCandidateLikVec.pop_back();
        deleteLikObject(lik_to_del);
    }
    omp_unset_lock(&mutex);

}

/**********************************************************************************************/
void ChromosomeNumberOptimizer::optimize(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds){
    int maxNumOfModels;
    SingleProcessPhyloLikelihood* lik;
    SingleProcessPhyloLikelihood* minAICcLik;
    
    std::map<int, vector<pair<uint, int>>> bestModelSharedParams;
    fixedParams_ = fixedParams;
    sharedParams_ = ChromEvolOptions::sharedParameters_;
    bestModelSharedParams = sharedParams_;
    optimizeBaseNumber_ = false;
    // if we should optimize it for at least one model, we will set it to true. Otherwise it is fixed for all the models.
    for (uint i = 1; i <= static_cast<uint>(ChromEvolOptions::numOfModels_); i++){
        if (!(std::count(fixedParams_[i].begin(), fixedParams_[i].end(), ChromosomeSubstitutionModel::BASENUM))){
            optimizeBaseNumber_ = true;
        }
    }


    ((ChromEvolOptions::maxNumOfModels_ == 1) && (ChromEvolOptions::heterogeneousModel_)) ? (maxNumOfModels = (static_cast<int>((tree_->getAllLeavesNames()).size())-1)) : (maxNumOfModels = ChromEvolOptions::maxNumOfModels_);
    vector<uint> candidateShiftNodesIds;
    if (ChromEvolOptions::heterogeneousModel_){
        getValidCandidatesForShift(candidateShiftNodesIds, ChromEvolOptions::minCladeSize_, ChromEvolOptions::numOfModels_);

    }
    
    vector <unsigned int> baseNumCandidates;
    getBaseNumCandidates(baseNumCandidates, baseNumberUpperBound_);
    uint numOfShifts = ChromEvolOptions::numOfModels_;

    std::map<int, std::vector<std::pair<uint, int>>> sharedParams = ChromEvolOptions::sharedParameters_;
    vectorOfLikelohoods_.reserve(numOfPoints);
    //vectorOfContexts_.reserve(numberOfModels);

    if (seed != 0){
        RandomTools::setSeed(static_cast<long>(seed));
    }
    bool deltaAICcImproved = true;
    bool firstIteration = true;

    while((deltaAICcImproved) && (numOfShifts <= (size_t)maxNumOfModels)){
        if ((candidateShiftNodesIds.size() == 0) && (maxNumOfModels > 1)){
            break;
        }
        
        bool improvedModelFound = false;
        if (firstIteration){
            initLikelihoods(modelParams, parsimonyBound, rateChange, numOfPoints, fixedRootFreqPath, fixedParams, mapModelNodesIds, numOfShifts, &sharedParams);
            optimizeMultiProcessModel(&sharedParams, &fixedParams, numOfPoints_, numOfIterations_, baseNumberUpperBound_, 0, 0);
            // leave only the best one
            clearVectorOfLikelihoods(1);
            //minAICcLik = vectorOfLikelohoods_[0];
            firstIteration = false;
            if ((prevModelsAICcLikValues_.empty()) && (ChromEvolOptions::heterogeneousModel_)){
                size_t numOfFixedParams = getNumberOfFixedParams(vectorOfLikelohoods_[0], fixedParams_); 
                double AICc = calculateModelSelectionCriterion(vectorOfLikelohoods_[0], numOfFixedParams);
                std::pair<double, double> AICcAndLik(AICc, vectorOfLikelohoods_[0]->getValue());
                prevModelsAICcLikValues_[ChromEvolOptions::numOfModels_] = AICcAndLik;
                getParameterNamesAndValues(vectorOfLikelohoods_[0], ChromEvolOptions::numOfModels_);
                prevModelsRootFrequencies_[ChromEvolOptions::numOfModels_] = getRootFrequencies(vectorOfLikelohoods_[0]);
            }
            numOfShifts ++;
            continue;

        }
        
        lik = vectorOfLikelohoods_.back();
        minAICcLik = lik;
        auto prevLik = minAICcLik;

        //vectorOfLikelohoods_.pop_back();     
        size_t numOfFixedParams = getNumberOfFixedParams(minAICcLik, fixedParams); 
        double initialAICc = calculateModelSelectionCriterion(minAICcLik, numOfFixedParams);
        uint minDetaAICcNode;
        double minAICc = initialAICc;
        std::cout << "*** *** *** Starting considering " << numOfShifts << " shifts *** *** ***" << std::endl;
        for (size_t i = 0; i < candidateShiftNodesIds.size(); i++){
            std::cout << "Candidate shifting node: N" << candidateShiftNodesIds[i] << " ..." << std::endl;
            uint prevShift = static_cast<uint>(lik->getSubstitutionProcess().getModelNumberForNode(candidateShiftNodesIds[i]));
            updateSharedParameters(sharedParams, prevShift, numOfShifts);
            fixedParams[numOfShifts] = fixedParams[prevShift];
            vectorOfLikelohoods_.pop_back();
            getNewLikObject(lik, candidateShiftNodesIds[i], &sharedParams_, &sharedParams, numOfPoints, fixedParams, parsimonyBound);
            optimizeMultiProcessModel(&sharedParams, &fixedParams, numOfPointsNextRounds_, numOfIterationsNextRounds_, baseNumberUpperBound_, 0, 0);
            clearVectorOfLikelihoods(1);
            auto candidateLik = vectorOfLikelohoods_[0];
            //optimizeModelParameters(candidateLik, ChromEvolOptions::tolerance_, ChromEvolOptions::maxIterations_, baseNumCandidates, &sharedParams, &fixedParams);
            numOfFixedParams = getNumberOfFixedParams(candidateLik, fixedParams); 
            double AICc_candidate =  calculateModelSelectionCriterion(candidateLik, numOfFixedParams);
            SingleProcessPhyloLikelihood* likToDel;
            // DEBUG //
            std::cout << ChromEvolOptions::modelSelectionCriterion_ <<" of the model before the shift:  " << initialAICc << std::endl;
            std::cout << ChromEvolOptions::modelSelectionCriterion_ << " of the candidate model : " << AICc_candidate << std::endl;
            if ((initialAICc - AICc_candidate > ChromEvolOptions::deltaAICcThreshold_) && (AICc_candidate < minAICc)){
                minDetaAICcNode = candidateShiftNodesIds[i];
                likToDel = minAICcLik;
                minAICcLik = candidateLik;
                minAICc = AICc_candidate;
                bestModelSharedParams = sharedParams;
                improvedModelFound = true;
                if (likToDel != prevLik){
                    deleteLikObject(likToDel);

                }
                
                std::cout << "** Best model is updated **" << std::endl;
            }else{
                likToDel = candidateLik;
                clearVectorOfLikelihoods(0);
                vectorOfLikelohoods_.push_back(minAICcLik);
                std::cout << "** Best model remains **" << std::endl;
            }
            printLikParameters(vectorOfLikelohoods_[0], 1, 0);
            
            //deleteLikObject(likToDel);
        }
        if (improvedModelFound){
            deltaAICcImproved = true;
            sharedParams_ = bestModelSharedParams;
            fixedParams_ = fixedParams;
            prevModelsPartitions_[numOfShifts] = prevModelsPartitions_[numOfShifts-1];
            prevModelsPartitions_[numOfShifts].push_back(minDetaAICcNode);
            std::pair<double, double> AICcAndLik(minAICc, minAICcLik->getValue());
            prevModelsAICcLikValues_[numOfShifts] = AICcAndLik;
            prevModelsRootFrequencies_[numOfShifts] = getRootFrequencies(minAICcLik);
            getParameterNamesAndValues(minAICcLik, numOfShifts);
            numOfShifts ++;
        }else{
            deltaAICcImproved = false;
           
        }   
        if (vectorOfLikelohoods_[0]->getSubstitutionProcess().getNumberOfModels() > 1){
            
            if (improvedModelFound){
                candidateShiftNodesIds.clear();
                getValidCandidatesForShift(candidateShiftNodesIds, ChromEvolOptions::minCladeSize_, numOfShifts-1);


            }
            
            //candidateShiftNodesIds.erase(std::remove(candidateShiftNodesIds.begin(), candidateShiftNodesIds.end(), minDetaAICcNode), candidateShiftNodesIds.end());

        }
        if (prevLik != minAICcLik){
            deleteLikObject(prevLik);
        }
        // sanity check //
       
        
    }
     std::cout << "*** Final best model: " << vectorOfLikelohoods_[0]->getValue() << std::endl;


}
/**********************************************************************************************/
size_t ChromosomeNumberOptimizer::getValidCandidatesForShiftRec(uint nodeId, std::vector<uint> &candidateShiftNodesIds, int minCladeSize, vector<uint> &shifting_nodes, bool shifting_node){
    auto sons = tree_->getSons(tree_->getNode(nodeId)); 
    size_t size = 0; 
    for (size_t i = 0;  i < sons.size(); i++){       
        uint sonIndex = tree_->getNodeIndex(sons[i]);
        if (tree_->isLeaf(sons[i])){
            size += 1;
        }else if (std::count(shifting_nodes.begin(), shifting_nodes.end(), sonIndex)){
            size += 1;
        }else{
            size+= getValidCandidatesForShiftRec(sonIndex, candidateShiftNodesIds, minCladeSize, shifting_nodes, false);

        }

    }
    if (((int)size >= minCladeSize) && (!shifting_node)){
        candidateShiftNodesIds.push_back(nodeId);
    }
    return size;
    

}
/**********************************************************************************************/
void ChromosomeNumberOptimizer::getValidCandidatesForShift(std::vector<uint> &candidateShiftNodesIds, int minCladeSize, uint numOfShifts){
    vector<shared_ptr<PhyloNode>> nodes = tree_->getAllNodes();
    auto &shifting_nodes = prevModelsPartitions_[numOfShifts];
    
    for (size_t j = 0; j < shifting_nodes.size(); j++){
        getValidCandidatesForShiftRec(shifting_nodes[j], candidateShiftNodesIds, minCladeSize, shifting_nodes, true);

    }

    
    // for (size_t i = 0; i < nodes.size(); i++){
    //     if (tree_->isLeaf(nodes[i])){
    //         continue;
    //     }
    //     if (tree_->getRootIndex() == tree_->getNodeIndex(nodes[i])){
    //         continue;
    //     }
    //     if (std::find(ChromEvolOptions::initialModelNodes_.begin(), ChromEvolOptions::initialModelNodes_.end(), tree_->getNodeIndex(nodes[i])) != ChromEvolOptions::initialModelNodes_.end()){
    //         continue;
    //     }
    //     auto leavesUnderNode = tree_->getLeavesUnderNode(nodes[i]);
    //     if (leavesUnderNode.size() >= (size_t)minCladeSize){
    //         candidateShiftNodesIds.push_back(tree_->getNodeIndex(nodes[i]));
    //     }
    // }

}

/*************************************************************************************/
uint ChromosomeNumberOptimizer::getNumberOfParametersPerParamType(int paramType, vector<int> &funcTypes){
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
/*********************************************************
* Functions for the Backward phase
**********************************************************/

void ChromosomeNumberOptimizer::getMapOfMergedModels(std::map<uint, uint> &mapOfModels, std::map<uint, vector<uint>> &modelsToMerge, SingleProcessPhyloLikelihood* prevLik){
    auto prevNumOfModels = static_cast<uint>(prevLik->getSubstitutionProcess().getNumberOfModels());
    vector<uint> modelsAfterMerge;
    vector<uint> removedModels;
    auto it = modelsToMerge.begin();
    while(it != modelsToMerge.end()){
        auto modelsPerClusterToRemove = modelsToMerge[it->first];
        for (size_t i = 0; i < modelsPerClusterToRemove.size(); i++){
            removedModels.push_back(modelsPerClusterToRemove[i]);
        }
        it++;
    }
    for (uint i = 1; i <= prevNumOfModels; i++){
        if (std::find(removedModels.begin(), removedModels.end(), i) == removedModels.end()){
            modelsAfterMerge.push_back(i);
        }
    }
    
    for (uint i = 1; i <= modelsAfterMerge.size(); i++){
        mapOfModels[modelsAfterMerge[i-1]] = i;
    }

}

void ChromosomeNumberOptimizer::mergeModels(std::map<uint, vector<uint>> modelsToMerge, SingleProcessPhyloLikelihood* lik, std::map<uint, vector<int>> &fixedParams, std::map<int, std::vector<std::pair<uint, int>>> &updatedSharedParams, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParamsPrevModel, std::map<uint, uint> &modelNums, std::map<uint, uint> &baseNumberBounds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams){
    getMapOfMergedModels(modelNums, modelsToMerge, lik);
    // updating shared and fixed parameters
    updateSharedParametersBackwards(sharedParams_, updatedSharedParams, modelNums);
    auto itModelMap = modelNums.begin();
    while (itModelMap != modelNums.end()){
        fixedParams[modelNums[itModelMap->first]] = fixedParams_[itModelMap->first];
        itModelMap ++;
    }

    // setting nodes partitions
    
    getMapOfModelAndNodeIdsBackward(mapModelNodesIds, lik, modelNums, modelsToMerge);

    // getting model parameters
    
    
    auto itModel = modelNums.begin();
    while (itModel != modelNums.end()){
        modelParams[modelNums[itModel->first]] = modelParamsPrevModel[itModel->first];
        auto branchProcess = lik->getSubstitutionProcess().getModel(itModel->first);
        baseNumberBounds[modelNums[static_cast<uint>(itModel->first)]] = std::dynamic_pointer_cast<const ChromosomeSubstitutionModel>(branchProcess)->getMaxChrRange();
        itModel++;
    }

}
////////////////////////////////////////////////////////////////////////////////////////
void ChromosomeNumberOptimizer::updateSharedParametersBackwards(std::map<int, std::vector<std::pair<uint, int>>> &sharedParams, std::map<int, std::vector<std::pair<uint, int>>> &updatedSharedParams, std::map<uint, uint> &mapOfModels){

    auto paramNumIt = sharedParams.begin();
    // removing the models that have to be merged
    while(paramNumIt != sharedParams.end()){
        int paramNum = paramNumIt->first;
        auto modelsAndTypes = sharedParams[paramNum];
        for (size_t i = 0; i < modelsAndTypes.size(); i++){
            uint model = sharedParams[paramNum][i].first;
            int type = sharedParams[paramNum][i].second;
            if (mapOfModels.find(model) != mapOfModels.end()){            
                std::pair<uint, int> modelAndType(mapOfModels[model], type);
                updatedSharedParams[paramNum].push_back(modelAndType);
            }
        }
        paramNumIt++;

    }
    return;

}

void ChromosomeNumberOptimizer::getMapOfModelAndNodeIdsBackward(std::map<uint, vector<uint>> &mapModelNodesIds, SingleProcessPhyloLikelihood* lik, std::map<uint, uint> &modelsMap, std::map<uint, vector<uint>> &modelsToMerge, uint rootId){
    uint numOfModels = static_cast<uint>(lik->getSubstitutionProcess().getNumberOfModels());
    if (rootId){
        mapModelNodesIds[modelsMap[1]].push_back(rootId);
    }
    // setting a map which will map the models that should be merged to the model they are merged to
    std::map<uint, uint> mapMergedToRemained;
    auto itMerged = modelsToMerge.begin();
    while(itMerged != modelsToMerge.end()){
        auto removedModels = modelsToMerge[itMerged->first];
        for (size_t j = 0; j < removedModels.size(); j++){
            auto removedModel = removedModels[j];
            mapMergedToRemained[removedModel] = itMerged->first;
        }      
        itMerged ++;
    }
    for (uint i = 1; i <= numOfModels; i++){
        auto vectorOfNodes = lik->getSubstitutionProcess().getNodesWithModel(i);
        auto it = modelsMap.find(i);
        uint model;
        if (it == modelsMap.end()){
            model = modelsMap[mapMergedToRemained[i]];
        }else{
            model = modelsMap[i];
        }
        for (size_t j = 0; j < vectorOfNodes.size(); j++){
            mapModelNodesIds[model].push_back(vectorOfNodes[j]);
        }

    }
}

pair<int, std::map<int, std::vector<double>>> ChromosomeNumberOptimizer::getMeanParameters(std::map<uint, pair<int, std::map<int, std::vector<double>>>> &prevModelParams, vector<uint> modelsToMerge, uint modelForMerge){
    std::pair<int, std::map<int, std::vector<double>>> meanValues;
    vector<int> baseNumbers;
    for (size_t i = 0; i < modelsToMerge.size(); i++){
        auto modelParams = prevModelParams[modelsToMerge[i]];
        if (i == 0){
            baseNumbers.push_back(modelParams.first);
            meanValues.second = modelParams.second;
        }else{
            baseNumbers.push_back(modelParams.first);
            auto it = modelParams.second.begin();
            while(it != modelParams.second.end()){
                for (size_t j = 0; j < modelParams.second[it->first].size(); j++){
                    meanValues.second[it->first][j] += modelParams.second[it->first][j];
                }               
                it ++;
            }

        }

    }
    auto modelParamsFirstModel = prevModelParams[modelForMerge];
    baseNumbers.push_back(modelParamsFirstModel.first);
    auto itTypeFirstModel = modelParamsFirstModel.second.begin();
    while (itTypeFirstModel != modelParamsFirstModel.second.end()){
        for (size_t j = 0; j < modelParamsFirstModel.second[itTypeFirstModel->first].size(); j++){
            meanValues.second[itTypeFirstModel->first][j] += modelParamsFirstModel.second[itTypeFirstModel->first][j];
        }
        itTypeFirstModel++;
    }
    sort(baseNumbers.begin(), baseNumbers.end());
    size_t midIndex = static_cast<size_t>(floor((double)(baseNumbers.size())/2));
    meanValues.first = baseNumbers[midIndex];
    auto itMeanVals = meanValues.second.begin();
    while(itMeanVals != meanValues.second.end()){
        for (size_t j = 0; j < meanValues.second[itMeanVals->first].size(); j++){
            meanValues.second[itMeanVals->first][j] /= (double)(modelsToMerge.size()+1);

        }
       
        itMeanVals ++;
    }
    return meanValues;

}

void ChromosomeNumberOptimizer::optimizeBackwards(double maxParsimony, bool parallel){
    backwardPhaseStarted_ = true;
    SingleProcessPhyloLikelihood* finalLikForward = vectorOfLikelohoods_[0];
    std::map<int, vector<std::pair<uint, int>>> sharedParams = sharedParams_;
    std::map<uint, vector<int>> fixedParams = fixedParams_;

    uint numOfClusters = static_cast<uint>(finalLikForward->getSubstitutionProcess().getNumberOfModels());
    if (numOfClusters <= 2){
        return;
    }
    auto finalLikBackward = finalLikForward;
    numOfShiftsForward_ = static_cast<uint>(finalLikBackward->getSubstitutionProcess().getNumberOfModels());
    size_t numOfFixedParams = getNumberOfFixedParams(finalLikForward, fixedParams_); 
    double AICc_best = calculateModelSelectionCriterion(finalLikForward, numOfFixedParams);
    bool areThereAreStillClusters = (numOfClusters > 0);
    while(areThereAreStillClusters){
        finalLikBackward = vectorOfLikelohoods_[0];
        std::map<std::pair<uint, uint>, double> pairsOfLikelihoods;
        uint numOfModels = static_cast<uint>(finalLikBackward->getSubstitutionProcess().getNumberOfModels());

        std::vector<std::pair<uint, uint>> pairsOfModels;
        // TODO: should also check sometimes combinations with the backward model
        for (uint i = 1; i <= numOfModels-1; i++){
            for (uint j = i+1; j <= numOfModels; j++){
                std::pair<uint, uint> pairOfModels(i, j);
                if (!isModelADirectSubtreeOfAnother(finalLikBackward, i, j)){
                    pairsOfModels.push_back(pairOfModels);

                }                
            }
        }
        time_t t1;
        time(&t1);
        time_t t2;
        if (parallel){
            omp_lock_t mutex;
            omp_init_lock(&mutex);
            #pragma omp parallel for schedule(dynamic)
            for (size_t i = 0; i < pairsOfModels.size(); i++){
                std::map<uint, vector<uint>> modelsToBeMerged;
                modelsToBeMerged[pairsOfModels[i].first].push_back(pairsOfModels[i].second);
                optimizeMergedModels(finalLikBackward, modelsToBeMerged, &pairsOfLikelihoods, maxParsimony, &mutex);

            }
            omp_destroy_lock(&mutex);

        }else{
            for (size_t i = 0; i < pairsOfModels.size(); i++){
                std::map<uint, vector<uint>> modelsToBeMerged;
                modelsToBeMerged[pairsOfModels[i].first].push_back(pairsOfModels[i].second);
                optimizeMergedModels(finalLikBackward, modelsToBeMerged, &pairsOfLikelihoods, maxParsimony, 0);

            }

        }
        time(&t2);
        std::cout <<"**** **** running time of the backward procedure is: "<< (t2-t1) <<endl;

        time_t t3;
        time(&t3);
        time_t t4;


        UndirectedGraph* G = new UndirectedGraph();
        for (size_t i = 1; i <= finalLikBackward->getSubstitutionProcess().getNumberOfModels(); i++){
            Vertex* modelNode = new Vertex(static_cast<uint>(i));
            G->addNewOrphanVertex(modelNode);

        }
        auto itEdges = pairsOfLikelihoods.begin();
        while(itEdges != pairsOfLikelihoods.end()){
            auto AICc_candidate = pairsOfLikelihoods[itEdges->first];
            if (AICc_best - AICc_candidate > ChromEvolOptions::deltaAICcThreshold_){
                G->addEdgeBetweenTwoNodes(itEdges->first.first, itEdges->first.second,  AICc_candidate);

            }
            
            itEdges ++;
        }

        std::map<uint, vector<Vertex*>> clusters;
        std::map<uint, vector<std::pair<uint, uint>>> maxEdges;
        G->countNumOfVerticesInClusters(clusters, maxEdges);
        auto fullyConectedClusters = G->findFullyConnectedClusters(clusters);
        areThereAreStillClusters = G->hasEdges();
        if (!areThereAreStillClusters){
            break;
        }
        // now checking the fully connected clusters
        // Iterate over each cluster, and merge the relevant models
        auto clusterRoots = G->getClusterRoots();
        std::map<uint, std::vector<uint>> rootAndVerticesToMerge;
        //SingleProcessPhyloLikelihood* better_lik = 0;
        for (size_t i = 0; i < clusterRoots.size(); i++){
            vector<SingleProcessPhyloLikelihood*> perClusterOfModelsLikVec;
            auto fullyConnected = fullyConectedClusters[clusterRoots[i]];
            //std::vector<uint> clusterNodes = G->getNodesIds(clusterRoots[i]);
            std::vector<uint> modelsToBeMerged;
            if (fullyConnected){
                for (size_t j = 0; j < clusters[clusterRoots[i]].size(); j++){
                    if (clusters[clusterRoots[i]][j]->getId() != clusterRoots[i]){
                        modelsToBeMerged.push_back(clusters[clusterRoots[i]][j]->getId());
                    }
                    
                }
                rootAndVerticesToMerge[clusterRoots[i]] = modelsToBeMerged;
                               
            }else{
                // bestEdge is a vector of best edges with the best AICc score
                // for now I will just use the first one (a random choice)        
                auto bestEdge = maxEdges[clusterRoots[i]];
                uint firstModel;
                uint secondModel;
                if (bestEdge[0].first < bestEdge[0].second){
                    firstModel = bestEdge[0].first;
                    secondModel = bestEdge[0].second;

                }else{
                    firstModel = bestEdge[0].second;
                    secondModel = bestEdge[0].first;
                }
                modelsToBeMerged.push_back(secondModel);
                rootAndVerticesToMerge[firstModel] = modelsToBeMerged;
            }


        }
        auto likToDel = finalLikBackward;
        vectorOfLikelohoods_.pop_back();
        mergeMultipleModelClusters(finalLikBackward, rootAndVerticesToMerge, maxParsimony);
        numOfFixedParams = getNumberOfFixedParams(vectorOfLikelohoods_[0], fixedParams_); 
        AICc_best = calculateModelSelectionCriterion(vectorOfLikelohoods_[0], numOfFixedParams);
        deleteLikObject(likToDel);
        delete G;
        time(&t4);
        std::cout <<"**** **** running time of the graph construction is: "<< (t4-t3) <<endl;
        // Delete all the unnecessary likelihood objects from the pairs map
        // HERE: TODO!! Is it relevant???

        

    }

}
void ChromosomeNumberOptimizer::optimizeMergedModels(SingleProcessPhyloLikelihood* prevLik, std::map<uint, vector<uint>> &modelsToBeMerged, std::map<std::pair<uint, uint>, double>* pairsOfLikelihoods, double parsimonyBound, omp_lock_t* mutex){
    vector<SingleProcessPhyloLikelihood*> perPairOfModelsLikVec;
    std::pair<uint, uint> pairOfMergedModels;
    uint numOfMerged = 0;

    auto itMergedModels = modelsToBeMerged.begin();
    while (itMergedModels != modelsToBeMerged.end()){
        if (pairsOfLikelihoods){
            pairOfMergedModels = std::pair<uint, uint>(itMergedModels->first, modelsToBeMerged[itMergedModels->first][0]);
        }
        numOfMerged += static_cast<uint>(modelsToBeMerged[itMergedModels->first].size());
        sort(modelsToBeMerged[itMergedModels->first].begin(), modelsToBeMerged[itMergedModels->first].end());
        itMergedModels++;
    }
    uint numOfModels = static_cast<uint>(prevLik->getSubstitutionProcess().getNumberOfModels())-numOfMerged;
    std::map<uint, vector<int>> fixedParams;
    std::map<int, std::vector<std::pair<uint, int>>> updatedSharedParams;
    std::map<uint, std::vector<uint>> mapModelNodesIds;
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> modelParams;
    std::map<uint, uint> mapOfModels;
    std::map<uint, uint> baseNumberBounds;
    std::string textToPrint = "";

    uint numOfModelsPrevModel = static_cast<uint>(prevLik->getSubstitutionProcess().getNumberOfModels());
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, 0, prevLik, &sharedParams_);   
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> modelParamsPrevModel = getMapOfParamsForComplexModel(prevLik, typeWithParamNames, numOfModelsPrevModel);
    if (mutex){
        omp_set_lock(mutex);
    }
    mergeModels(modelsToBeMerged, prevLik, fixedParams, updatedSharedParams, mapModelNodesIds, modelParamsPrevModel, mapOfModels, baseNumberBounds, modelParams);
    if (mutex){
        omp_unset_lock(mutex);

    }
    optimizeFirstRound(&updatedSharedParams, fixedParams, parsimonyBound, mapModelNodesIds, modelParams, numOfModels, numOfPointsNextRounds_, numOfIterationsNextRounds_, perPairOfModelsLikVec, &textToPrint, &baseNumberBounds, &mapOfModels, &modelParamsPrevModel, &modelsToBeMerged, mutex);
    optimizeMultiProcessModel(&updatedSharedParams, &fixedParams, numOfPointsNextRounds_, numOfIterationsNextRounds_, baseNumberBounds, &perPairOfModelsLikVec, &textToPrint, mutex);

    size_t numOfFixedParams = getNumberOfFixedParams(perPairOfModelsLikVec[0], fixedParams); 
    double AICc = calculateModelSelectionCriterion(perPairOfModelsLikVec[0], numOfFixedParams);
    auto modelToDel = perPairOfModelsLikVec.back();
    if (mutex){
        omp_set_lock(mutex);

    }
    if (pairsOfLikelihoods){
        (*pairsOfLikelihoods)[pairOfMergedModels] = AICc;

    }
    std::cout << "\tMerging models: " << pairOfMergedModels.first << ", " << pairOfMergedModels.second << std::endl;
    std::cout << textToPrint << std::endl; 
    //std::cout << "optimized log likelihood is: " <<  perPairOfModelsLikVec[0]->getValue() << std::endl;
    printLikParameters(perPairOfModelsLikVec[0], 1, 0);
    std::cout << "Final " << ChromEvolOptions::modelSelectionCriterion_ << " is: " << AICc << std::endl;
    deleteLikObject(modelToDel);
    if (mutex){
        omp_unset_lock(mutex);

    }
    

}
//************************************************************************************/
void ChromosomeNumberOptimizer::mergeMultipleModelClusters(SingleProcessPhyloLikelihood* prevLik, std::map<uint, vector<uint>> &rootAndVerticesToMerge, double maxParsimony){


    std::map<uint, vector<int>> fixedParams;
    std::map<int, std::vector<std::pair<uint, int>>> updatedSharedParams;
    std::map<uint, std::vector<uint>> mapModelNodesIds;
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> modelParams;
    std::map<uint, uint> mapOfModels;
    std::map<uint, uint> baseNumberBounds;

    uint numOfModelsPrevModel = static_cast<uint>(prevLik->getSubstitutionProcess().getNumberOfModels());
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, 0, prevLik, &sharedParams_);   
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> modelParamsPrevModel = getMapOfParamsForComplexModel(prevLik, typeWithParamNames, numOfModelsPrevModel);

    mergeModels(rootAndVerticesToMerge, prevLik, fixedParams, updatedSharedParams, mapModelNodesIds, modelParamsPrevModel, mapOfModels, baseNumberBounds, modelParams);
    auto numOfModels = static_cast<uint>(mapOfModels.size());
    optimizeFirstRound(&updatedSharedParams, fixedParams, maxParsimony, mapModelNodesIds, modelParams, numOfModels, numOfPointsNextRounds_, numOfIterationsNextRounds_, vectorOfLikelohoods_, 0, &baseNumberBounds, &mapOfModels, &modelParamsPrevModel, &rootAndVerticesToMerge, 0);
    optimizeMultiProcessModel(&updatedSharedParams, &fixedParams, numOfPointsNextRounds_, numOfIterationsNextRounds_, baseNumberBounds, 0, 0, 0);
    sharedParams_ = updatedSharedParams;
    fixedParams_ = fixedParams;
    // get here the AICc?

}
/*********************************************************************************************/
void ChromosomeNumberOptimizer::fillWithFathers(vector<uint> &fathers, vector<uint> &nodes){
    for (size_t i = 0; i < nodes.size(); i++){
        if (tree_->getRootIndex() != nodes[i]){
            auto fatherNode = tree_->getFatherOfNode(tree_->getNode(nodes[i]));
            uint father = tree_->getNodeIndex(fatherNode);
            fathers.push_back(father);
        }
    }
}
/*********************************************************************************************/
bool ChromosomeNumberOptimizer::isModelADirectSubtreeOfAnother(SingleProcessPhyloLikelihood* lik, size_t indexModel1, size_t indexModel2){
    auto vectorOfNodesM1 = lik->getSubstitutionProcess().getNodesWithModel(indexModel1);
    auto vectorOfNodesM2 = lik->getSubstitutionProcess().getNodesWithModel(indexModel2);
    if (indexModel1 == 1){
        vectorOfNodesM1.push_back(tree_->getRootIndex());
    }
    vector<uint> fathersM1;
    vector<uint> fathersM2;
    fillWithFathers(fathersM1, vectorOfNodesM1);
    set<uint> setFathersM1(fathersM1.begin(), fathersM1.end());
    set<uint> setNodesM1(vectorOfNodesM1.begin(),vectorOfNodesM1.end());
    set<uint> setNodesM2(vectorOfNodesM2.begin(), vectorOfNodesM2.end());
    set<uint> intersectWithM1;
    set<uint> intersectWithM2;
    set_intersection(setFathersM1.begin(), setFathersM1.end(), setNodesM2.begin(), setNodesM2.end(),
                 std::inserter(intersectWithM2, intersectWithM2.begin()));
    set_intersection(setFathersM1.begin(), setFathersM1.end(), setNodesM1.begin(), setNodesM1.end(),
                 std::inserter(intersectWithM1, intersectWithM1.begin()));
    // If M1 is contained within M2, the fathers should be either contained 
    // within M1 or M2             
    if ((intersectWithM2.size() > 0) && ((intersectWithM1.size() + intersectWithM2.size()) == setFathersM1.size())){
        return true;
    }
    intersectWithM1.clear();
    intersectWithM2.clear();
    // a symmetrical check for fathers of M2 to check whether M2 is contained within M1.
    fillWithFathers(fathersM2, vectorOfNodesM2);
    set<uint> setFathersM2(fathersM2.begin(), fathersM2.end());
    set_intersection(setFathersM2.begin(), setFathersM2.end(), setNodesM2.begin(), setNodesM2.end(),
                 std::inserter(intersectWithM2, intersectWithM2.begin()));
    set_intersection(setFathersM2.begin(), setFathersM2.end(), setNodesM1.begin(), setNodesM1.end(),
                 std::inserter(intersectWithM1, intersectWithM1.begin()));
    if ((intersectWithM1.size() > 0) && ((intersectWithM1.size() + intersectWithM2.size()) == setFathersM2.size())){
        return true;
    }
    return false;

}