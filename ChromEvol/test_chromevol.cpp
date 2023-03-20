//standard libraries
#include <string>
#include <vector>
#include <iostream>
#include <time.h>
#include <unistd.h>
#include <cstdlib>
//from bpp-core
#include <Bpp/Version.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/BppApplication.h>




//from bpp-seq
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>




#include "ChromEvolOptions.h"
#include "ChromosomeSubstitutionModel.h"
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/NewLikelihood/MarginalAncestralReconstruction.h>
#include <Bpp/Phyl/NewLikelihood/JointMLAncestralReconstruction.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowNumeric.h>
#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/RateAcrossSitesSubstitutionProcess.h>
#include "ChrFasta.h"
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>




using namespace bpp;
using namespace std;


void setMaxChrNum(unsigned int maxNumberOfChr, int &maxChrNum){
    if (maxChrNum < 0){
        maxChrNum = maxNumberOfChr + std::abs(maxChrNum);
    }else{
        if ((int)maxNumberOfChr > maxChrNum){
            maxChrNum = maxNumberOfChr;
        }
    }

}
/****************************************************************************/
void setMinChrNum(unsigned int minNumberOfChr, int &minChrNum){
    if (minChrNum < 0){
        if (minNumberOfChr == 1){
            minChrNum = 0;
            std::cout << "Warning !!!! minChrNum_ should be at least 1!!" << std::endl;
            std::cout << "The mininal chromosome number was determined to be 1" << std::endl;
        }
        minChrNum = minNumberOfChr - std::abs(minChrNum);
    }else{
        if ((int)minNumberOfChr < minChrNum){
            minChrNum = minNumberOfChr;
        }

    }
}
/*************************************************************************************************************/
VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc, ChromosomeAlphabet* alphaInitial, int &minChrNum, int &maxChrNum, ChromosomeAlphabet** alphabet){
    size_t numOfSequences = vsc->getNumberOfSequences();
    vector <string> sequenceNames = vsc->getSequencesNames();
    *alphabet = new ChromosomeAlphabet(minChrNum,maxChrNum);
        // fill with composite values
    if (alphaInitial->getNumberOfCompositeStates() > 0){
        const std::map <int, std::map<int, double>> compositeStates = alphaInitial->getCompositeStatesMap();
        std::map <int, std::map<int, double>>::const_iterator it = compositeStates.begin();
        while (it != compositeStates.end()){
            int compositeState = it->first;
            std::string charComposite = alphaInitial->intToChar(compositeState);
            (*alphabet)->setCompositeState(charComposite);
            it++;
        }
    }
    VectorSiteContainer* resized_alphabet_site_container = new VectorSiteContainer(*alphabet);
    for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = vsc->getSequence(sequenceNames[i]);
        BasicSequence new_seq = BasicSequence(seq.getName(), seq.getChar(0), *alphabet);
        resized_alphabet_site_container->addSequence(new_seq);

    }
    return resized_alphabet_site_container;
}

/****************************************************************************/
VectorSiteContainer* getCharacterData (const string& path, uint &numberOfUniqueStates, std::map<uint, uint> &chrRange, std::map<uint, int> &baseNum, int numOfModels, ChromosomeAlphabet** alphabet){
    ChromosomeAlphabet* alphaInitial = new ChromosomeAlphabet(1, 500);
    VectorSequenceContainer* initialSetOfSequences = ChrFasta::readSequencesFromFile(path, alphaInitial);
    size_t numOfSequences = initialSetOfSequences->getNumberOfSequences();
    vector <string> sequenceNames = initialSetOfSequences->getSequencesNames();

    unsigned int maxNumberOfChr = 1; //the minimal number of chromosomes cannot be zero
    unsigned int minNumOfChr = 500;

    std::vector <int> UniqueCharacterStates;
    cout<<"vector size is "<< UniqueCharacterStates.size()<<endl;
    for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = initialSetOfSequences->getSequence(sequenceNames[i]);
        int character = seq.getValue(0);
        if (character == -1){
            continue;
        }
        if (character == 501){
            continue;
        }
        // if it is a composite state
        if (character > 501){
            const std::vector<int> compositeCharacters = alphaInitial->getSetOfStatesForAComposite(character);
            for (size_t j = 0; j < compositeCharacters.size(); j++){
                if ((unsigned int) compositeCharacters[j] > maxNumberOfChr){
                    maxNumberOfChr = compositeCharacters[j];
                }
                if ((unsigned int) compositeCharacters[j] < minNumOfChr){
                    minNumOfChr = compositeCharacters[j];
                }
                
            }
            continue;
        }

        if (!std::count(UniqueCharacterStates.begin(), UniqueCharacterStates.end(), character)){
            UniqueCharacterStates.push_back(character);
        }
        if ((unsigned int) character > maxNumberOfChr){
            maxNumberOfChr = character;
        }
        if ((unsigned int) character < minNumOfChr){
            minNumOfChr = character;
        }

    }
    numberOfUniqueStates = (unsigned int)UniqueCharacterStates.size() + alphaInitial->getNumberOfCompositeStates();
    uint chrRangeNum = maxNumberOfChr - minNumOfChr;
    for (uint j = 1; j <= static_cast<uint>(numOfModels); j++){      
        if (baseNum[j] != IgnoreParam){
            if (baseNum[j] > (int)chrRangeNum){
                chrRange[j] = baseNum[j] + 1;
            }else{
                chrRange[j] = chrRangeNum;
            }
        }
    }
    cout <<"Number of unique states is " << numberOfUniqueStates <<endl;
    int maxChrNum = -10;
    int minChrNum = -1;

    setMaxChrNum(maxNumberOfChr, maxChrNum);
    setMinChrNum(minNumOfChr, minChrNum);

    VectorSiteContainer* vsc = resizeAlphabetForSequenceContainer(initialSetOfSequences, alphaInitial, minChrNum, maxChrNum, alphabet);
    delete initialSetOfSequences;
    delete alphaInitial;
    return vsc;
}


/****************************************************************************/
void rescale_tree(PhyloTree* tree, double chrRange){
    bool rooted = tree->isRooted();
    if (!rooted){
        throw Exception("The given input tree is unrooted. Tree must be rooted!\n");
    }

    double treeLength = tree->getTotalLength();
    double scale_tree_factor = chrRange/treeLength;
    tree->scaleTree(scale_tree_factor);

}
/****************************************************************************/
PhyloTree* getTree(const string& path, double treeLength){
    Newick reader;
    PhyloTree* tree = reader.readPTree(path);
    rescale_tree(tree, treeLength);
    return tree;

}
/****************************************************************************/
void deleteLikObject(SingleProcessPhyloLikelihood* lik_to_del){
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
/****************************************************************************/
SingleProcessPhyloLikelihood* createLikInstance(VectorSiteContainer* vsc, PhyloTree* tree, uint &numberOfUniqueStates, std::map<uint, uint> &chrRange, std::map<uint, int> &baseNum, int baseNumber, int numOfModels, std::map<uint, std::vector<uint>> &mapModelNodesIds, ChromosomeAlphabet* alphabet){


    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    ParametrizablePhyloTree* parTree = new ParametrizablePhyloTree(*tree);
    std::map<int, std::vector<double>> mapOfParamValues;
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::GAIN)].push_back(2);
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::LOSS)].push_back(2);
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DUPL)].push_back(3);
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)].push_back(1.3);
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)].push_back(0.1);
    std::vector<int> rateChangeType;
    for (size_t i = 0; i < 5; i++){
        rateChangeType.push_back(static_cast<int>(ChromosomeNumberDependencyFunction::CONSTANT));

    }


    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, mapOfParamValues, baseNumber, chrRange[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType);
    subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree);


    // adding models
    for (uint i = 1; i <= numOfModels; i++){
        if (i > 1){
            chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, mapOfParamValues,baseNumber, chrRange[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType);
        }   
        subProSim->addModel(chrModel, mapModelNodesIds[i]);
    }

    SubstitutionProcess* nsubPro= subProSim->clone();
    Context* context = new Context();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc->clone(), *nsubPro, true);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;

}
/****************************************************************************/


int main() {
    time_t t1;
    time(&t1);
    time_t t2;
    time_t t3;
    time_t t4;
    const string characterFilePath = "/home/anat/Docs/Crinum/counts.fasta";
    const string treeFilePath = "/home/anat/Docs/Crinum/tree.newick";
    ChromosomeAlphabet* alphabet; 

    int baseNumber = 11;
    uint numberOfUniqueStates;
    std::map<uint, uint> chrRange;
    std::map<uint, int> baseNum;
    baseNum[1] = baseNumber;
    int numOfModels = 1;
    std::map<uint, std::vector<uint>> mapModelNodesIds;
    VectorSiteContainer* vsc = getCharacterData (characterFilePath, numberOfUniqueStates, chrRange, baseNum, numOfModels, &alphabet);
    PhyloTree* tree = getTree(treeFilePath, (double)numberOfUniqueStates);
    // adding nodes corresponding to the model
    auto nodes = tree->getAllNodes();
    for (size_t i = 0; i < nodes.size(); i++){
        uint nodeId = tree->getNodeIndex(nodes[i]);
        if (nodeId == tree->getRootIndex()){
            continue;
        }else{
            mapModelNodesIds[1].push_back(nodeId);
        }
    }
    sleep(60);
    time(&t2);
    std::cout <<"Total running time before initialization: "<< static_cast<int>(t2-t1) <<endl;



    std::vector<SingleProcessPhyloLikelihood*> likVector;
    size_t numberOfItems = 1000;
    for (size_t i = 0; i < numberOfItems; i++){
        auto lik = createLikInstance(vsc, tree, numberOfUniqueStates, chrRange, baseNum, baseNumber, numOfModels, mapModelNodesIds, alphabet);
        likVector.push_back(lik);
    }
    sleep(60);
    time(&t3);
    std::cout <<"Total running time after initialization: "<< static_cast<int>(t3-t1) <<endl;
    size_t numOfObjects = likVector.size();
    

    for (size_t i = 0; i < numOfObjects; i++){
        SingleProcessPhyloLikelihood* lik_to_del =likVector.back(); 
        likVector.pop_back();
        deleteLikObject(lik_to_del);

    }

    delete vsc;
    delete tree;
    delete alphabet;
    sleep(60);
    time(&t4);
    std::cout <<"Total running time after deletion: "<< static_cast<int>(t4-t1) <<endl;

    return 0;
}
