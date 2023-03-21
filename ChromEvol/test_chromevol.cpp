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
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>

#include "sys/types.h"
#include "sys/sysinfo.h"




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
SubstitutionProcess* createProcess(PhyloTree* tree, const NucleicAlphabet* alphabet,  std::map<uint, std::vector<uint>> &mapModelNodesIds){
    auto rootFreqs = std::make_shared<GCFrequencySet>(alphabet);
    auto model = std::make_shared<T92>(alphabet, 3., .1);
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(4, 1.0);
    ParametrizablePhyloTree* parTree = new ParametrizablePhyloTree(*tree);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree);
    subProSim->addModel(model, mapModelNodesIds[1]);
    SubstitutionProcess* nsubPro= subProSim->clone();
    return nsubPro;

}

/****************************************************************************/
SingleProcessPhyloLikelihood* createLikInstanceNonChromevol(VectorSiteContainer* vsc, PhyloTree* tree, const NucleicAlphabet* alphabet,  std::map<uint, std::vector<uint>> &mapModelNodesIds){
    
    auto rootFreqs = std::make_shared<GCFrequencySet>(alphabet);
    auto model = std::make_shared<T92>(alphabet, 3., .1);
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(4, 1.0);
    ParametrizablePhyloTree* parTree = new ParametrizablePhyloTree(*tree);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree);
    subProSim->addModel(model, mapModelNodesIds[1]);
    SubstitutionProcess* nsubPro= subProSim->clone();
    Context* context = new Context();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc->clone(), *nsubPro, false);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;
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
void test_memory(){
    std::vector<int*> vec;
    // 1000000000
    for (int i = 0; i < 10000000; i++){
        vec.push_back(new int(i));
    }
    system("free -h");
    size_t vec_size = vec.size();
    for (size_t i = 0; i < vec_size; i++){
        int* toDel = vec.back();
        delete toDel;
        vec.pop_back();

    }
    vec.clear();
    vec.shrink_to_fit();
    return;
}
/****************************************************************************/
void test_chromevol_memory(){
    time_t t1;
    time(&t1);
    time_t t2;
    time_t t3;
    time_t t4;
    const string characterFilePath = "/home/anat/Docs/Crinum/counts.fasta";
    const string treeFilePath = "/home/anat/Docs/Crinum/tree.newick";
    //system("free -h");
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
    time(&t2);
    std::cout <<"Total running time before initialization: "<< static_cast<int>(t2-t1) <<endl;
    std::vector<SingleProcessPhyloLikelihood*> likVector;
    // //std::vector<SubstitutionProcess*> likVector;
    // std::vector<VectorSiteContainer*> likVector;
    size_t numberOfItems = 10000;
    for (size_t i = 0; i < numberOfItems; i++){
        auto lik = createLikInstance(vsc, tree, numberOfUniqueStates, chrRange, baseNum, baseNumber, numOfModels, mapModelNodesIds, alphabet);
    //     //auto lik = createLikInstanceNonChromevol(vsc, tree, alphabet, mapModelNodesIds);
    //     //auto process = createProcess(tree, alphabet,  mapModelNodesIds);
    //     auto sites = vsc->clone();


        auto lik_val = lik->getValue();
        std::cout << "lik value is: " << lik_val << std::endl;
        likVector.push_back(lik);
    //     likVector.push_back(sites);
    //     //likVector.push_back(process);
    }
    sleep(5);
    time(&t3);
    std::cout <<"Total running time after initialization: "<< static_cast<int>(t3-t1) <<endl;
    size_t numOfObjects = likVector.size();
    system("free -h");

    

    for (size_t i = 0; i < numOfObjects; i++){
        SingleProcessPhyloLikelihood* lik_to_del =likVector.back(); 
        //SubstitutionProcess* processToDel = likVector.back(); 
        //VectorSiteContainer* vsc_to_del = likVector.back(); 
        likVector.pop_back();
        //delete vsc_to_del;
        //delete processToDel;
        deleteLikObject(lik_to_del);

    }


    delete vsc;
    delete tree;
    sleep(15);
    time(&t4);
    std::cout <<"Total running time after deletion: "<< static_cast<int>(t4-t1) <<endl;
    //system("free -h");

}
/****************************************************************************/
void test(){
    std::vector<int> vec;
    for (int i = 0; i < 10; i++){
        vec.push_back(i);
    }
    for (size_t i = 0; i < 5; i++){
        vec.pop_back();
    }
    vec.shrink_to_fit();
    for (size_t i = 0; i < vec.size(); i++){
        std::cout << vec[i] << std::endl;
    }
}
/****************************************************************************/


int main() {
    // Newick reader;
    // PhyloTree* tree = reader.parenthesisToPhyloTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);", false, "", false, false);
    // const NucleicAlphabet* alphabet = &AlphabetTools::DNA_ALPHABET;


    system("free -h");
    test_memory();
    //test();
    //test_chromevol_memory();

    //vec = vector<int*>();
    //delete vec;
    //vec.shrink_to_fit();


    // std::vector<int>* vec = new vector<int>();
    // for (int i = 0; i < 1000000000; i++){
    //     vec->push_back(i);
    // }
    // size_t vec_size = vec->size();
    // for (size_t i = 0; i < vec_size; i++){
    //     vec->pop_back();
    // }
    // delete vec;
    //sleep(10);
    //system("free -h");
    // time_t t1;
    // time(&t1);
    // time_t t2;
    // time_t t3;
    // time_t t4;
    // // const st ring characterFilePath = "/home/anat/Docs/Crinum/counts.fasta";
    // // const string treeFilePath = "/home/anat/Docs/Crinum/tree.newick";
    // // system("free -h");
    // // ChromosomeAlphabet* alphabet; 

    // // int baseNumber = 11;
    // // uint numberOfUniqueStates;
    // // std::map<uint, uint> chrRange;
    // // std::map<uint, int> baseNum;
    // // baseNum[1] = baseNumber;
    // // int numOfModels = 1;
    // std::map<uint, std::vector<uint>> mapModelNodesIds;
    // // VectorSiteContainer* vsc = getCharacterData (characterFilePath, numberOfUniqueStates, chrRange, baseNum, numOfModels, &alphabet);
    // // PhyloTree* tree = getTree(treeFilePath, (double)numberOfUniqueStates);
    // // adding nodes corresponding to the model
    // auto nodes = tree->getAllNodes();
    // for (size_t i = 0; i < nodes.size(); i++){
    //     uint nodeId = tree->getNodeIndex(nodes[i]);
    //     if (nodeId == tree->getRootIndex()){
    //         continue;
    //     }else{
    //         mapModelNodesIds[1].push_back(nodeId);
    //     }
    // }

    // VectorSiteContainer* vsc = new VectorSiteContainer(alphabet);
    // vsc->addSequence(BasicSequence("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", alphabet));
    // vsc->addSequence(BasicSequence("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", alphabet));
    // vsc->addSequence(BasicSequence("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", alphabet));
    // vsc->addSequence(BasicSequence("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", alphabet));

    // sleep(5);
    // time(&t2);
    // std::cout <<"Total running time before initialization: "<< static_cast<int>(t2-t1) <<endl;
    // system("free -h");



    // //std::vector<SingleProcessPhyloLikelihood*> likVector;
    // //std::vector<SubstitutionProcess*> likVector;
    // std::vector<VectorSiteContainer*> likVector;
    // size_t numberOfItems = 30000;
    // for (size_t i = 0; i < numberOfItems; i++){
    //     //auto lik = createLikInstance(vsc, tree, numberOfUniqueStates, chrRange, baseNum, baseNumber, numOfModels, mapModelNodesIds, alphabet);
    //     //auto lik = createLikInstanceNonChromevol(vsc, tree, alphabet, mapModelNodesIds);
    //     //auto process = createProcess(tree, alphabet,  mapModelNodesIds);
    //     auto sites = vsc->clone();


    //     //auto lik_val = lik->getValue();
    //     //std::cout << "lik value is: " << lik_val << std::endl;
    //     //likVector.push_back(lik);
    //     likVector.push_back(sites);
    //     //likVector.push_back(process);
    // }
    // sleep(5);
    // time(&t3);
    // std::cout <<"Total running time after initialization: "<< static_cast<int>(t3-t1) <<endl;
    // size_t numOfObjects = likVector.size();
    // system("free -h");

    

    // for (size_t i = 0; i < numOfObjects; i++){
    //     //SingleProcessPhyloLikelihood* lik_to_del =likVector.back(); 
    //     //SubstitutionProcess* processToDel = likVector.back(); 
    //     VectorSiteContainer* vsc_to_del = likVector.back(); 
    //     likVector.pop_back();
    //     delete vsc_to_del;
    //     //delete processToDel;
    //     //deleteLikObject(lik_to_del);

    // }


    // delete vsc;
    // delete tree;
    sleep(15);
    //std::cout <<"Total running time after deletion: "<< static_cast<int>(t4-t1) <<endl;
    system("free -h");

    return 0;
    
}
