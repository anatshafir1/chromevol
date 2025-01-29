#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/JointPhyloLikelihood.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/Character/CharacterSubstitutionModel.h>
#include <Bpp/Phyl/Model/Character/SingleRateModel.h>
#include <Bpp/Phyl/Model/Character/RatePerEntryModel.h>
#include <Bpp/Phyl/Model/Character/RatePerExitModel.h>
#include <Bpp/Phyl/Model/Character/RatePerPairSymModel.h>
#include <Bpp/Phyl/Model/Character/RatePerPairModel.h>
#include "ChromosomeNumberOptimizer.h"
#include "JointTraitChromosomeLikelihood.h"
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include "BrownianMotionAncestralReconstruction.h"
#include "BrownianMotionLikelihood.h"
#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>

using namespace bpp;
using namespace std;


// void printNodeRec(std::shared_ptr<PhyloTree> tree, uint nodeId, bool names){
//     if ((!(tree->isLeaf(nodeId)))){
//         auto sons = tree->getSons(nodeId);
//         if (names){
//             std::cout << "current node is " << (tree->getNode(nodeId))->getName();
//             std::cout << " id is " << nodeId << std::endl;

//         }else{
//             std::cout << "current node is " << nodeId << std::endl;
//         }
        
//         std::cout << "\tSons are: " << std::endl;
//         for (size_t i = 0; i < sons.size(); i++){
//             shared_ptr<PhyloBranch> branch=tree->getEdgeToFather(sons[i]);
//             if (names){
//                 std::cout << "\t\t" << (tree->getNode(sons[i]))->getName() << " branch length is " << branch->getLength();
//                 std::cout << " id is " << sons[i] << std::endl;

//             }else{
//                 std::cout << "\t\t" <<sons[i] << " branch length is " << branch->getLength() << std::endl;
//             }
            
//         }
//         for (size_t i = 0; i < sons.size(); i++){
//             printNodeRec(tree, sons[i], names);
//         }
//     }
// }
// /******************************************************/
// void printTree(std::shared_ptr<PhyloTree> tree, bool names){
//     uint rootNodeId = tree->getRootIndex();
//     printNodeRec(tree, rootNodeId, names);

// }
// /******************************************************/
// std::shared_ptr<PhyloTree> createExpectedMappingHistory(std::shared_ptr<LikelihoodCalculationSingleProcess> lik, uint mappingsNum, StochasticMapping* stm){
  
//     //ChromEvolOptions::NumOfSimulations_);
//     stm->generateStochasticMapping();
//     std::vector<std::shared_ptr<PhyloTree>> mappings;
//     for (size_t i = 0; i < mappingsNum; i++){
//         auto mappingTree = stm->createMappingHistoryTree(i);
//         std::cout << "Mapping tree " << i << std::endl;
//         mappings.push_back(mappingTree);
//         printTree(mappingTree, true);
//         std::cout << "***|***|***|***|****" << std::endl;
//     }
//     std::cout << "Expected mapping tree: " << std::endl;
        
//     auto expected_mapping = stm->generateExpectedMapping(mappings);
//     printTree(expected_mapping, true);
//     return expected_mapping;
// }

// std::map <uint, std::vector<uint>> getNodesForEachModel(std::shared_ptr<PhyloTree> expectedMapping, StochasticMapping* stm){
//   std::map <uint, std::vector<uint>> nodeModels;
//   auto nodes = expectedMapping->getAllNodes();
//   for (size_t i = 0; i < nodes.size(); i ++){
//     uint nodeId = expectedMapping->getNodeIndex(nodes[i]);
//     int nodeState = stm->getNodeState(nodes[i]);
//     if (!(expectedMapping->isLeaf(nodeId))){
//       auto sons = expectedMapping->getSons(nodeId);
//       for (size_t j = 0; j < sons.size(); j++){
//         // add +1 so that it will work for the general function that constructs the heterogeneous model
//         nodeModels[static_cast<uint>(nodeState)+1].push_back(sons[j]);
//       }
//     }
//   }
//   return nodeModels;
// }


void getTraitNodes(std::shared_ptr<PhyloTree> tree, vector<uint> &modelNodes){
  auto nodes = tree->getAllNodes();
  for (size_t i = 0; i < nodes.size(); i++){
    auto nodeId = tree->getNodeIndex(nodes[i]);
    if (nodeId == tree->getRootIndex()){
      continue;
    }
    modelNodes.push_back(nodeId);
  }
}
// std::shared_ptr<TwoParameterBinarySubstitutionModel> initBinaryTraitModel(double pi0, double mu, std::shared_ptr<PhyloTree> phyloTree){
//   std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*phyloTree);
//   const BinaryAlphabet* alphabet = new BinaryAlphabet();
//   std::shared_ptr<DiscreteDistribution> rdistTrait = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
//   auto twoParamModel = std::make_shared<TwoParameterBinarySubstitutionModel>(alphabet,mu,pi0);
//   RowMatrix<double> pijt = twoParamModel->getPij_t(0.3);
//   MatrixTools::printForR(pijt,"Pijt",std::cout);
//   auto& Q = twoParamModel->getGenerator();
//   MatrixTools::printForR(Q,"Q",std::cout);
//   RowMatrix<double> pijt_dt = twoParamModel->getdPij_dt(0.3);
//   MatrixTools::printForR(pijt_dt,"Pijt_dt",std::cout);
//   RowMatrix<double> d2Pij_dt2 = twoParamModel->getd2Pij_dt2(0.3);
//   MatrixTools::printForR(d2Pij_dt2,"d2Pij_dt2",std::cout);
//   std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);
//   std::vector<uint> modelNodesTrait;
//   getTraitNodes(phyloTree, modelNodesTrait);
//   std::shared_ptr<VectorSiteContainer> sites = std::make_shared<VectorSiteContainer>(alphabet);
// 	sites->addSequence(BasicSequence("S1", "1", alphabet));
// 	sites->addSequence(BasicSequence("S2", "1", alphabet));
// 	sites->addSequence(BasicSequence("S3", "0", alphabet));
// 	sites->addSequence(BasicSequence("S4", "0", alphabet));
// 	sites->addSequence(BasicSequence("S5", "1", alphabet));
  
      
      
//   subProT->addModel(std::shared_ptr<TwoParameterBinarySubstitutionModel>(twoParamModel->clone()), modelNodesTrait);
//   Context context;
//   auto likT = std::make_shared<LikelihoodCalculationSingleProcess>(context, *sites->clone(), *(subProT->clone()));
//   auto phyloT = new SingleProcessPhyloLikelihood(context, likT);

//   std::cout << "Initial likelihood is: " << phyloT->getValue() << std::endl;
//   std::cout << "Model parameters are: " << std::endl;
//   auto substitutionModelParams = phyloT->getSubstitutionModelParameters();
//   for (size_t i = 0; i < substitutionModelParams.size(); i++){
//     std::cout << "\t" << substitutionModelParams[i].getName() << " value is " << substitutionModelParams[i].getValue() << std::endl;
//   }
//   return twoParamModel;

// }

// std::shared_ptr<CharacterSubstitutionModel> initIntegerTraitLikelihood(vector<double> freqVals, double globalrate, std::shared_ptr<PhyloTree> phyloTree){
//   std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*phyloTree);
//   const IntegerAlphabet* traitAlpha = new IntegerAlphabet(1);
//   std::shared_ptr<DiscreteDistribution> rdistTrait = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));

//   shared_ptr<IntegerFrequencySet> freqs = make_shared<FullIntegerFrequencySet>(traitAlpha, freqVals);
//   std::shared_ptr<CharacterSubstitutionModel> characterModel = std::make_shared<SingleRateModel>(traitAlpha, freqs, false);
  
//   characterModel->setFrequencySet(*freqs);
//   characterModel->setParameterValue("global_rate", globalrate);
//   RowMatrix<double> pijt = characterModel->getPij_t(0.3);
//   MatrixTools::printForR(pijt,"Pijt",std::cout);
//   auto& Q = characterModel->getGenerator();
//   MatrixTools::printForR(Q,"Q",std::cout);
//   RowMatrix<double> pijt_dt = characterModel->getdPij_dt(0.3);
//   MatrixTools::printForR(pijt_dt,"Pijt_dt",std::cout);
//   RowMatrix<double> d2Pij_dt2 = characterModel->getd2Pij_dt2(0.3);
//   MatrixTools::printForR(d2Pij_dt2,"d2Pij_dt2",std::cout);
//   std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);
//   std::vector<uint> modelNodesTrait;
//   getTraitNodes(phyloTree, modelNodesTrait);
//   std::shared_ptr<VectorSiteContainer> sites = std::make_shared<VectorSiteContainer>(traitAlpha);
// 	sites->addSequence(BasicSequence("S1", "1", traitAlpha));
// 	sites->addSequence(BasicSequence("S2", "1", traitAlpha));
// 	sites->addSequence(BasicSequence("S3", "0", traitAlpha));
// 	sites->addSequence(BasicSequence("S4", "0", traitAlpha));
// 	sites->addSequence(BasicSequence("S5", "1", traitAlpha));
      
//   subProT->addModel(std::shared_ptr<CharacterSubstitutionModel>(characterModel->clone()), modelNodesTrait);
//   Context context;
//   auto likT = std::make_shared<LikelihoodCalculationSingleProcess>(context, *sites->clone(), *(subProT->clone()));
//   auto phyloT = new SingleProcessPhyloLikelihood(context, likT);

//   std::cout << "Initial likelihood is: " << phyloT->getValue() << std::endl;
//   std::cout << "Model parameters are: " << std::endl;
//   auto substitutionModelParams = phyloT->getSubstitutionModelParameters();
//   for (size_t i = 0; i < substitutionModelParams.size(); i++){
//     std::cout << "\t" << substitutionModelParams[i].getName() << " value is " << substitutionModelParams[i].getValue() << std::endl;
//   }
//   return characterModel;

// }

// void testBinary(shared_ptr<PhyloTree> pTree){
//    //std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*pTree);
//   double globalRate = 2.0;
//   double pi0 = 0.25;
//   vector<double> freqs;
//   freqs.push_back(pi0);
//   freqs.push_back(1-pi0);
//   std::cout << "******" << std::endl;
//   std::cout << "Binary Model:" << std::endl;
//   std::shared_ptr<TwoParameterBinarySubstitutionModel> binary_model = initBinaryTraitModel(pi0, globalRate, pTree);
//   std::cout << "******" << std::endl;
//   std::shared_ptr<CharacterSubstitutionModel> integer_model = initIntegerTraitLikelihood(freqs, globalRate, pTree);

//   for (size_t i = 0; i < 2; i++){
//     for (size_t j = 0; j < 2; j++){
//       if (abs(binary_model->Pij_t(i, j, 0.3) - integer_model->Pij_t(i, j, 0.3)) > 0.00001){
//         std::cout << "ERROR!!!" << std::endl;
//         return;
//       }
//       if (abs(binary_model->dPij_dt(i, j, 0.3) - integer_model->dPij_dt(i, j, 0.3)) > 0.00001){
//         std::cout << "ERROR!!!" << std::endl;
//         return;
//       }
//       if (abs(binary_model->d2Pij_dt2(i, j, 0.3) - integer_model->d2Pij_dt2(i, j, 0.3)) > 0.00001){
//         std::cout << "ERROR!!!" << std::endl;
//         return;
//       }
//     }
//   }
//   std::cout << "Success" << std::endl; 

// }
// /****************************************************************************************/
void setParametersToNewTraitModel(std::map<string, double> &traitModelParams, std::shared_ptr<CharacterSubstitutionModel> traitModel, shared_ptr<IntegerFrequencySet> freqs){
  traitModel->setFrequencySet(*freqs);
  auto it = traitModelParams.begin();
  while (it != traitModelParams.end()){
    auto paramName = it->first;
    string prefix = "pi";
    if (paramName.compare(0, prefix.length(), prefix) != 0){
      traitModel->setParameterValue(paramName, traitModelParams[paramName]);
      
    }
    it ++;
  }
}
// /****************************************************************************************/
// void getBaseNumCandidates(std::vector <unsigned int> &baseNumCandidates, std::shared_ptr<VectorSiteContainer> seqData){
//     size_t numOfSequences = seqData->getNumberOfSequences();
//     unsigned int minRange = 0;
//     vector <string> sequenceNames = seqData->getSequenceNames();
//     for (size_t i = 0; i < numOfSequences; i++){
//         if (i == numOfSequences-1){
//             continue;
//         }
//         BasicSequence seq1 = seqData->getSequence(sequenceNames[i]);
//         int chrNum1 = seq1.getValue(0);
//         if (chrNum1 == -1){
//             continue;
//         }
//         for (size_t j = i + 1; j < numOfSequences; j++){
//             BasicSequence seq2 = seqData->getSequence(sequenceNames[j]);
//             int chrNum2 = seq2.getValue(0);
//             if (chrNum2 == -1){
//                 continue;
//             }
//             unsigned int chrRange = (unsigned int)(std::abs(chrNum1 - chrNum2));
//             if (chrRange < lowerLimitBaseNumber){
//                 continue;
//             }

//             if (!std::count(baseNumCandidates.begin(), baseNumCandidates.end(), chrRange)){
//                 if (minRange == 0){
//                     minRange = chrRange;
//                 }else{
//                     if (chrRange < minRange){
//                         minRange = chrRange;
//                     }
//                 }
//                 baseNumCandidates.push_back(chrRange);

//             }

//         }
//     }
//     if (minRange > lowerLimitBaseNumber){
//         for (unsigned int i = lowerLimitBaseNumber; i < minRange; i++){
//             baseNumCandidates.push_back(i);
//         }

//     }

// }
// /****************************************************************************************/
std::shared_ptr<CharacterSubstitutionModel> setTraitModel(const IntegerAlphabet* traitAlpha, shared_ptr<IntegerFrequencySet> freqset, std::string modelName){
    if (modelName == "singleRate"){
        return std::make_shared<SingleRateModel>(traitAlpha, freqset, false);
    }else if (modelName == "symPairRate"){
        return std::make_shared<RatePerPairSymModel>(traitAlpha, freqset, false);
    }else if (modelName == "pairRateModel"){
        return std::make_shared<RatePerPairModel>(traitAlpha, freqset, false);
    }else if (modelName == "ratePerEntry"){
        return std::make_shared<RatePerEntryModel>(traitAlpha, freqset, false);
    }else if (modelName == "ratePerExit"){
        return std::make_shared<RatePerExitModel>(traitAlpha, freqset, false);
    }else{
        throw Exception("ERROR!!! LikelihoodUtils::setTraitModel(): No such trait model exists!!!");
    }
}

// /****************************************************************************************/
// void testMulti1(std::map<string, double> traitModelParams, std::shared_ptr<PhyloTree> tree){
//   std::string modelName = "singleRate";
//   std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree.get());
//   std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT;
//   std::shared_ptr<NonHomogeneousSubstitutionProcess> nsubProT;
//   std::shared_ptr<LikelihoodCalculationSingleProcess> likT;
//   //bool traitWeightedFrequencies = true;
//   std::shared_ptr<BranchModel> traitModel;
//   std::vector<uint> modelNodesTrait;
//   std::shared_ptr<DiscreteDistribution> rdistTrait = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));  
//   Context contextT = Context();
//   const IntegerAlphabet* traitAlpha = new IntegerAlphabet(2);
//   vector<double> freqVals;
//   freqVals.push_back(0.3);
//   freqVals.push_back(0.6);
//   freqVals.push_back(0.1);
//   shared_ptr<IntegerFrequencySet> freqs = make_shared<FullIntegerFrequencySet>(traitAlpha, freqVals);
//   std::shared_ptr<CharacterSubstitutionModel> model = setTraitModel(traitAlpha, freqs, modelName);
  
//   //ratePerPairModel.setParameterValue("rate_1_0", rateVals(1, 0));
//   setParametersToNewTraitModel(traitModelParams, model, freqs);

//   traitModel = static_pointer_cast<BranchModel>(model);
//   subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);
//   getTraitNodes(tree, modelNodesTrait);
    
//   subProT->addModel(std::shared_ptr<CharacterSubstitutionModel>(model->clone()), modelNodesTrait);
//   nsubProT = std::shared_ptr<NonHomogeneousSubstitutionProcess>(subProT->clone());

  
//   Context context;
//   std::shared_ptr<VectorSiteContainer> vscTrait = std::make_shared<VectorSiteContainer>(traitAlpha);
// 	vscTrait->addSequence(BasicSequence("S1", "1", traitAlpha));
// 	vscTrait->addSequence(BasicSequence("S2", "2", traitAlpha));
// 	vscTrait->addSequence(BasicSequence("S3", "0", traitAlpha));
// 	vscTrait->addSequence(BasicSequence("S4", "0", traitAlpha));
// 	vscTrait->addSequence(BasicSequence("S5", "1", traitAlpha));
//   likT = std::make_shared<LikelihoodCalculationSingleProcess>(contextT, *vscTrait->clone(), *nsubProT);
//   auto phyloT = std::make_shared<SingleProcessPhyloLikelihood>(contextT, likT);
//   phyloT->getValue(); // calculate now because it will be needed for stochastic mapping

//   // get expected tree for the chromosome model
//   size_t numOfStochasticMappings = 3;
//   StochasticMapping* stm = new StochasticMapping(likT, numOfStochasticMappings);
//   std::map<uint, std::vector<size_t>> mlAncestors = JointPhyloLikelihood::getMLAncestralReconstruction(phyloT.get(), phyloT->getParameters(), 2);
//   stm->setMLAncestors(&mlAncestors);
//   auto treeChr = createExpectedMappingHistory(likT, numOfStochasticMappings, stm);
//   //auto treeChr = stm->createExpectedMappingHistory(numOfStochasticMappings);
//   // creating the chromosome number model
//   auto rdistC = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
//   std::map<uint, vector<uint>> mapModelNodesIds = getNodesForEachModel(treeChr, stm);
//   delete stm;
//   std::shared_ptr<ParametrizablePhyloTree> parTreeChr;
//   std::vector<std::shared_ptr<ChromosomeSubstitutionModel>> models;
//   SubstitutionProcess* subProcess;
//   const ChromosomeAlphabet* alphabetChr = new ChromosomeAlphabet(3,20);
//   auto rdist2 = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
//   auto parTree2 =  std::make_shared<ParametrizablePhyloTree>(*treeChr);
//   std::map<int, std::vector<double>> mapOfParamValues;
//   mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::GAIN)].push_back(2);
//   mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::LOSS)].push_back(2);
//   mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DUPL)].push_back(3);
//   mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)].push_back(1.3);
//   mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)].push_back(0.1);
//   uint chrRange = 15;
//   int baseNumber = 6;
//   std::vector<int> rateChangeType;
//   for (size_t i = 0; i < 5; i++){
//     rateChangeType.push_back(static_cast<int>(ChromosomeNumberDependencyFunction::CONSTANT));
//   }
//   std::shared_ptr<ChromosomeSubstitutionModel> chrModel1 = std::make_shared<ChromosomeSubstitutionModel>(alphabetChr, mapOfParamValues, baseNumber, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType);
//   std::shared_ptr<ChromosomeSubstitutionModel> chrModel2 = std::make_shared<ChromosomeSubstitutionModel>(alphabetChr, mapOfParamValues, baseNumber, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType);
//   std::shared_ptr<ChromosomeSubstitutionModel> chrModel3 = std::make_shared<ChromosomeSubstitutionModel>(alphabetChr, mapOfParamValues, baseNumber, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType);

//   std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim2;
//   subProSim2 = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdist2->clone()), parTree2);
//   auto nodesMap = getNodesForEachModel(treeChr, stm);
//   subProSim2->addModel(chrModel1, nodesMap[1]);
//   subProSim2->addModel(chrModel2, nodesMap[2]);
//   subProSim2->addModel(chrModel2, nodesMap[3]);
//   subProcess= subProSim2->clone();

//   auto modelColl=std::make_shared<SubstitutionProcessCollection>();
//   models.push_back(chrModel1);
//   models.push_back(chrModel2);
//   models.push_back(chrModel3);
  
//   modelColl->addModel(traitModel, 1);
//   for (size_t i = 0; i < models.size(); i++){
//     modelColl->addModel(models[i], i+2);

//   }


//   std::shared_ptr<DiscreteDistribution> rdistChr = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
  
//   modelColl->addDistribution(rdistTrait, 1);
//   modelColl->addDistribution(rdistChr, 2);

//   modelColl->addTree(parTree, 1);
//   modelColl->addTree(parTree2, 2);
//   vector<uint> vP1m1 = modelNodesTrait;

//   map<size_t, Vuint> mModBr1;
//   mModBr1[1]=vP1m1;
 

//   modelColl->addSubstitutionProcess(1, mModBr1, 1, 1);
                                   
//   map<size_t, Vuint> mModBr2;

//   for (size_t i = 0; i < models.size(); i++){
//     mModBr2[i+2]= mapModelNodesIds[i+1];
//   }


//   modelColl->addSubstitutionProcess(2, mModBr2, 2, 2);

//   // Now setting the joint likelihood object
//     // Likelihoods
//   SubstitutionProcess* sP1c=nsubProT->clone();
//   SubstitutionProcess* sP2c=subProcess->clone(); // no need to clone again- it was already cloned.
//   delete subProcess;
//   auto pc(std::make_shared<PhyloLikelihoodContainer>(context, *modelColl));
//   std::shared_ptr<VectorSiteContainer> vscChr = std::make_shared<VectorSiteContainer>(alphabetChr);
// 	vscChr->addSequence(BasicSequence("S1-1", "3", alphabetChr));
// 	vscChr->addSequence(BasicSequence("S2-1", "6", alphabetChr));
// 	vscChr->addSequence(BasicSequence("S3-0", "18", alphabetChr));
// 	vscChr->addSequence(BasicSequence("S4-0", "12", alphabetChr));
// 	vscChr->addSequence(BasicSequence("S5-1", "6", alphabetChr));



//   auto lik1_j = std::make_shared<LikelihoodCalculationSingleProcess>(context, *(vscTrait->clone()), *sP1c);

//   pc->addPhyloLikelihood(1, new SingleProcessPhyloLikelihood(context, lik1_j));
    
//   auto lik2_j = std::make_shared<LikelihoodCalculationSingleProcess>(context, *(vscChr->clone()), *sP2c, true);

//   pc->addPhyloLikelihood(2, new SingleProcessPhyloLikelihood(context, lik2_j));
//   std::vector <unsigned int> baseNumberCandidates;
//   getBaseNumCandidates(baseNumberCandidates, vscChr);
//   bool ml= true;
//   //   JointTraitChromosomeLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool expectedHistory, bool weightedFrequencies, size_t numOfMappings, std::string traitModel, std::vector<int> &rateChangeType, VectorSiteContainer* chromosomeVsc, std::vector<unsigned int> &baseNumberCandidates, bool inCollection = true);
//   JointTraitChromosomeLikelihood* jl = new JointTraitChromosomeLikelihood(context, pc, true, true, numOfStochasticMappings, modelName, rateChangeType, vscChr->clone(), baseNumberCandidates, ml);
//   jl->setStochasticMappingTree(treeChr);
//   //jl->setSharedParams(sharedParams_);
//   //jl->setFixedParams(fixedParams_);
//   std::cout << "Initial likelihood is " << jl->getValue() << std::endl;
//   auto params = jl->getSubstitutionModelParameters();
//   std::cout <<"Parameters are:" << std::endl;
//   for (size_t i = 0; i < params.size(); i++){
//     std::cout <<"\t" << params[i].getName() << ": " << params[i].getValue() << std::endl;

//   }

// }
SubstitutionProcess* setChromosomeSubstitutionModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, std::map<int, vector<std::pair<uint,int>>>* updatedSharedParams, bool weightedRootFreqs, vector<std::shared_ptr<ChromosomeSubstitutionModel>>* models, std::shared_ptr<ParametrizablePhyloTree> parTree, vector<int> &rateChangeType){
    std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    if (parTree == nullptr){
        parTree = std::make_shared<ParametrizablePhyloTree>(*tree);

    }
    //std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
    std::map<uint, std::map<int, vector<string>>> mapOfParamsNamesPerModelType;
    LikelihoodUtils::setParamsNameInForMultiProcess(mapOfParamsNamesPerModelType, modelParams);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, modelParams[1].second, modelParams[1].first, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType, false);
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
            chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, modelParams[i].second, modelParams[i].first, baseNumberUpperBound[i], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType, false);
        }   
        subProSim->addModel(chrModel, mapModelNodesIds[i]);
        if (models){
            models->push_back(chrModel);
        }
    }
    //LikelihoodUtils::aliasParametersInSubstitutionProcess(mapOfParamsNamesPerModelType, updatedSharedParams, subProSim);


    SubstitutionProcess* nsubPro= subProSim->clone();
    auto nonhomoProcess = dynamic_cast<NonHomogeneousSubstitutionProcess*>(nsubPro);
    LikelihoodUtils::aliasParametersInSubstitutionProcess(mapOfParamsNamesPerModelType, updatedSharedParams, nonhomoProcess);
    return nsubPro;

}
std::map <uint, std::vector<uint>> getNodesForEachModel(std::shared_ptr<PhyloTree> expectedMapping, StochasticMapping* stm){
  std::map <uint, std::vector<uint>> nodeModels;
  auto nodes = expectedMapping->getAllNodes();
  for (size_t i = 0; i < nodes.size(); i ++){
    uint nodeId = expectedMapping->getNodeIndex(nodes[i]);
    int nodeState = stm->getNodeState(nodes[i]);
    if (!(expectedMapping->isLeaf(nodeId))){
      auto sons = expectedMapping->getSons(nodeId);
      for (size_t j = 0; j < sons.size(); j++){
        // add +1 so that it will work for the general function that constructs the heterogeneous model
        nodeModels[static_cast<uint>(nodeState)+1].push_back(sons[j]);
      }
    }
  }
  return nodeModels;
}


/**************************************************************************************************************************************/
uint getMaxBaseNumAmongModels(std::map<uint, uint> baseNumberUpperBound) {
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
/******************************************************************************************/
void printLikelihoodParameters(SingleProcessPhyloLikelihood* lik, unsigned int optimized, vector<string> paramsNames){

 
    for (int i = 0; i < (int)(paramsNames.size()); i++){
        if (paramsNames[i].find("Chromosome.baseNum_") != std::string::npos){
            std::cout << paramsNames[i] << " = " << (int)(lik->getParameter(paramsNames[i]).getValue()) << std::endl;
        }else{
            std::cout << paramsNames[i] << " = " << lik->getParameter(paramsNames[i]).getValue() << std::endl;

        }    
    }
    std::cout <<  "***" << std::endl;

}
/******************************************************************************************/
void optimizeTraitModel(SingleProcessPhyloLikelihood* lik, double tol, uint numOfIterations, std::vector<string> &traitParamNames, std::vector<string> &fixedParamsTrait){
    DerivableSecondOrder* f = lik;
    BrentOneDimension* optimizer = new BrentOneDimension(f);
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);
    // setting bracketing for Brent optimization
    optimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);
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

        }
        printLikelihoodParameters(lik, 1, traitParamNames);
        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
    }
    delete optimizer;

}

// /***************************************************************************************/
void getAllPossibleChrRanges(std::vector <unsigned int> &baseNumCandidates, std::shared_ptr<VectorSiteContainer> seqData) {
    // size_t numOfSequences = seqData->getNumberOfSequences();
    // unsigned int minRange = 0;
    // vector <string> sequenceNames = seqData->getSequenceNames();
    // for (size_t i = 0; i < numOfSequences; i++){
    //     if (i == numOfSequences-1){
    //         continue;
    //     }
    //     BasicSequence seq1 = seqData->getSequence(sequenceNames[i]);
    //     int chrNum1 = seq1.getValue(0);
    //     if (chrNum1 == -1){
    //         continue;
    //     }
    //     for (size_t j = i + 1; j < numOfSequences; j++){
    //         BasicSequence seq2 = seqData->getSequence(sequenceNames[j]);
    //         int chrNum2 = seq2.getValue(0);
    //         if (chrNum2 == -1){
    //             continue;
    //         }
    //         unsigned int chrRange = (unsigned int)(std::abs(chrNum1 - chrNum2));
    //         if (chrRange < 3){
    //             continue;
    //         }

    //         if (!std::count(baseNumCandidates.begin(), baseNumCandidates.end(), chrRange)){
    //             if (minRange == 0){
    //                 minRange = chrRange;
    //             }else{
    //                 if (chrRange < minRange){
    //                     minRange = chrRange;
    //                 }
    //             }
    //             baseNumCandidates.push_back(chrRange);
                

    //         }

    //     }
    // }
    // if (minRange > 3){
    //     for (unsigned int i = 3; i < minRange; i++){
    //         baseNumCandidates.push_back(i);
    //     }

    // }
    baseNumCandidates.push_back(3);
    baseNumCandidates.push_back(15);
    baseNumCandidates.push_back(9);
    baseNumCandidates.push_back(12);
    baseNumCandidates.push_back(6);

}
/******************************************************/
void getBaseNumCandidates(vector <unsigned int> &baseNumCandidates, std::map<uint, uint> &baseNumberUpperBounds, std::shared_ptr<VectorSiteContainer> seqData){
    //uint maxBaseNumCandidate = getMaxBaseNumAmongModels(baseNumberUpperBounds);
    getAllPossibleChrRanges(baseNumCandidates, seqData);

}
/******************************************************/
void testDefault(){
  double seedUb = 10000000;
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((S1:0.1,S2:0.1):0.3,S3:0.4):0.2,(S4:0.3,S5:0.3):0.3);", false, "", false, false));
  RandomTools::setSeed(static_cast<long int>(seedUb));

  std::map<string, double> traitModelParams;
  traitModelParams["global_rate"] =  2;
  std::string modelName = "singleRate";
  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT;
  std::shared_ptr<NonHomogeneousSubstitutionProcess> nsubProT;
  std::shared_ptr<LikelihoodCalculationSingleProcess> likT;

  Context contextT = Context();
  const IntegerAlphabet* traitAlpha = new IntegerAlphabet(1);
  vector<double> freqVals;
  freqVals.push_back(0.3);
  freqVals.push_back(0.7);
  


  std::shared_ptr<VectorSiteContainer> vscTrait = std::make_shared<VectorSiteContainer>(traitAlpha);
	vscTrait->addSequence(BasicSequence("S1", "1", traitAlpha));
	vscTrait->addSequence(BasicSequence("S2", "1", traitAlpha));
	vscTrait->addSequence(BasicSequence("S3", "0", traitAlpha));
	vscTrait->addSequence(BasicSequence("S4", "0", traitAlpha));
	vscTrait->addSequence(BasicSequence("S5", "1", traitAlpha));


  std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*pTree);

  //bool traitWeightedFrequencies = true;
  std::shared_ptr<BranchModel> traitModel;
  std::vector<uint> modelNodesTrait;
  std::shared_ptr<DiscreteDistribution> rdistTrait = make_shared<GammaDiscreteRateDistribution>(1, 1.0);

  shared_ptr<IntegerFrequencySet> freqs = make_shared<FullIntegerFrequencySet>(traitAlpha, freqVals);
  std::shared_ptr<CharacterSubstitutionModel> model = setTraitModel(traitAlpha, freqs, modelName);
  
  //ratePerPairModel.setParameterValue("rate_1_0", rateVals(1, 0));
  setParametersToNewTraitModel(traitModelParams, model, freqs);
  auto& Q = model->getGenerator();
  MatrixTools::printForR(Q,"Q",std::cout);

  traitModel = static_pointer_cast<BranchModel>(model);
  subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);
  getTraitNodes(pTree, modelNodesTrait);
    
  subProT->addModel(std::shared_ptr<CharacterSubstitutionModel>(model->clone()), modelNodesTrait);
  nsubProT = std::shared_ptr<NonHomogeneousSubstitutionProcess>(subProT->clone());
  string prefix  = traitModel->getNamespace();
  std::unordered_map<std::string, std::string> sharedTraitParams;
  int numOfConstraints = 0;
  //LikelihoodUtils::aliasTraitParams(std::shared_ptr<NonHomogeneousSubstitutionProcess> process, int &numOfTraitConstraints, string &prefix, std::unordered_map<string, string> &sharedTraitParams);
  LikelihoodUtils::aliasTraitParams(nsubProT, numOfConstraints, prefix, sharedTraitParams);
  auto traitData = vscTrait->clone();

  likT = std::make_shared<LikelihoodCalculationSingleProcess>(contextT, *traitData, *nsubProT, false);
  SingleProcessPhyloLikelihood* phyloT = new SingleProcessPhyloLikelihood(contextT, likT);
  std::cout << "likelihood is: " << phyloT->getValue() << std::endl;// calculate now because it will be needed for stochastic mapping
  
  auto params = phyloT->getSubstitutionModelParameters();
  std::cout <<"Parameters are:" << std::endl;
  for (size_t i = 0; i < params.size(); i++){
    std::cout <<"\t" << params[i].getName() << ": " << params[i].getValue() << std::endl;

  }
  // std::vector<string> parameterNames = params.getParameterNames();
  // vector<string> fixed_params_names;
  // optimizeTraitModel(phyloT, 0.1, 5, parameterNames, fixed_params_names);
  // std::cout << "final likelihood is: " << phyloT->getValue() << std::endl;
  delete traitAlpha;
  delete phyloT;


}
/******************************************************/
void testTraitWeighted(){
  double seedUb = 10000000;
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((S1:0.1,S2:0.1):0.3,S3:0.4):0.2,(S4:0.3,S5:0.3):0.3);", false, "", false, false));
  RandomTools::setSeed(static_cast<long int>(seedUb));

  std::map<string, double> traitModelParams;
  traitModelParams["global_rate"] = 2.0;
  std::string modelName = "singleRate";
  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT;
  std::shared_ptr<NonHomogeneousSubstitutionProcess> nsubProT;
  std::shared_ptr<LikelihoodCalculationSingleProcess> likT;

  Context contextT = Context();
  const IntegerAlphabet* traitAlpha = new IntegerAlphabet(1);
  vector<double> freqVals;
  freqVals.push_back(0.3);
  freqVals.push_back(0.7);
  


  std::shared_ptr<VectorSiteContainer> vscTrait = std::make_shared<VectorSiteContainer>(traitAlpha);
	vscTrait->addSequence(BasicSequence("S1", "1", traitAlpha));
	vscTrait->addSequence(BasicSequence("S2", "1", traitAlpha));
	vscTrait->addSequence(BasicSequence("S3", "0", traitAlpha));
	vscTrait->addSequence(BasicSequence("S4", "0", traitAlpha));
	vscTrait->addSequence(BasicSequence("S5", "1", traitAlpha));


  std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*pTree);

  //bool traitWeightedFrequencies = true;
  std::shared_ptr<BranchModel> traitModel;
  std::vector<uint> modelNodesTrait;
  std::shared_ptr<DiscreteDistribution> rdistTrait = make_shared<GammaDiscreteRateDistribution>(1, 1.0);

  shared_ptr<IntegerFrequencySet> freqs = make_shared<FullIntegerFrequencySet>(traitAlpha, freqVals);
  std::shared_ptr<CharacterSubstitutionModel> model = setTraitModel(traitAlpha, freqs, modelName);
  
  //ratePerPairModel.setParameterValue("rate_1_0", rateVals(1, 0));
  setParametersToNewTraitModel(traitModelParams, model, freqs);
  auto& Q = model->getGenerator();
  MatrixTools::printForR(Q,"Q",std::cout);

  traitModel = static_pointer_cast<BranchModel>(model);
  subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);
  getTraitNodes(pTree, modelNodesTrait);
    
  subProT->addModel(std::shared_ptr<CharacterSubstitutionModel>(model->clone()), modelNodesTrait);
  nsubProT = std::shared_ptr<NonHomogeneousSubstitutionProcess>(subProT->clone());
  string prefix  = traitModel->getNamespace();
  std::unordered_map<std::string, std::string> sharedTraitParams;
  int numOfConstraints = 0;
  //LikelihoodUtils::aliasTraitParams(std::shared_ptr<NonHomogeneousSubstitutionProcess> process, int &numOfTraitConstraints, string &prefix, std::unordered_map<string, string> &sharedTraitParams);
  LikelihoodUtils::aliasTraitParams(nsubProT, numOfConstraints, prefix, sharedTraitParams);
  auto traitData = vscTrait->clone();

  likT = std::make_shared<LikelihoodCalculationSingleProcess>(contextT, *traitData, *nsubProT, true);
  SingleProcessPhyloLikelihood* phyloT = new SingleProcessPhyloLikelihood(contextT, likT);
  std::cout << "likelihood is: " << phyloT->getValue() << std::endl;// calculate now because it will be needed for stochastic mapping
  
  auto params = phyloT->getSubstitutionModelParameters();
  std::cout <<"Parameters are:" << std::endl;
  for (size_t i = 0; i < params.size(); i++){
    std::cout <<"\t" << params[i].getName() << ": " << params[i].getValue() << std::endl;

  }
  std::vector<string> parameterNames = params.getParameterNames();
  vector<string> fixed_params_names;
  optimizeTraitModel(phyloT, 0.1, 5, parameterNames, fixed_params_names);
  std::cout << "final likelihood is: " << phyloT->getValue() << std::endl;
  delete traitAlpha;
  delete phyloT;


}
void testTraitFixed(){
  double seedUb = 10000000;
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((S1:0.1,S2:0.1):0.3,S3:0.4):0.2,(S4:0.3,S5:0.3):0.3);", false, "", false, false));
  RandomTools::setSeed(static_cast<long int>(seedUb));

  std::map<string, double> traitModelParams;
  traitModelParams["global_rate"] = 2.0;
  std::string modelName = "singleRate";
  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT;
  std::shared_ptr<NonHomogeneousSubstitutionProcess> nsubProT;
  std::shared_ptr<LikelihoodCalculationSingleProcess> likT;

  Context contextT = Context();
  const IntegerAlphabet* traitAlpha = new IntegerAlphabet(1);
  vector<double> freqVals;
  freqVals.push_back(0.3);
  freqVals.push_back(0.7);
  


  std::shared_ptr<VectorSiteContainer> vscTrait = std::make_shared<VectorSiteContainer>(traitAlpha);
	vscTrait->addSequence(BasicSequence("S1", "1", traitAlpha));
	vscTrait->addSequence(BasicSequence("S2", "1", traitAlpha));
	vscTrait->addSequence(BasicSequence("S3", "0", traitAlpha));
	vscTrait->addSequence(BasicSequence("S4", "0", traitAlpha));
	vscTrait->addSequence(BasicSequence("S5", "1", traitAlpha));


  std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*pTree);

  //bool traitWeightedFrequencies = true;
  std::shared_ptr<BranchModel> traitModel;
  std::vector<uint> modelNodesTrait;
  std::shared_ptr<DiscreteDistribution> rdistTrait = make_shared<GammaDiscreteRateDistribution>(1, 1.0);

  shared_ptr<IntegerFrequencySet> freqs = make_shared<FullIntegerFrequencySet>(traitAlpha, freqVals);
  std::shared_ptr<CharacterSubstitutionModel> model = setTraitModel(traitAlpha, freqs, modelName);
  
  //ratePerPairModel.setParameterValue("rate_1_0", rateVals(1, 0));
  setParametersToNewTraitModel(traitModelParams, model, freqs);
  auto& Q = model->getGenerator();
  MatrixTools::printForR(Q,"Q",std::cout);

  vector<double> rootFreqVals;
  rootFreqVals.push_back(1);
  rootFreqVals.push_back(0);
  std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(model->getStateMap(), false)), rootFreqVals);
  std::shared_ptr<FrequencySet> rootFrequenciesFixedRoot = std::shared_ptr<FrequencySet>(rootFreqsFixed->clone());
  
  traitModel = static_pointer_cast<BranchModel>(model);
  subProT = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree, rootFrequenciesFixedRoot);
  getTraitNodes(pTree, modelNodesTrait);
    
  subProT->addModel(std::shared_ptr<CharacterSubstitutionModel>(model->clone()), modelNodesTrait);
  nsubProT = std::shared_ptr<NonHomogeneousSubstitutionProcess>(subProT->clone());
  string prefix  = traitModel->getNamespace();
  std::unordered_map<std::string, std::string> sharedTraitParams;
  int numOfConstraints = 0;
  //LikelihoodUtils::aliasTraitParams(std::shared_ptr<NonHomogeneousSubstitutionProcess> process, int &numOfTraitConstraints, string &prefix, std::unordered_map<string, string> &sharedTraitParams);
  LikelihoodUtils::aliasTraitParams(nsubProT, numOfConstraints, prefix, sharedTraitParams);
  auto traitData = vscTrait->clone();

  likT = std::make_shared<LikelihoodCalculationSingleProcess>(contextT, *traitData, *nsubProT, false);
  SingleProcessPhyloLikelihood* phyloT = new SingleProcessPhyloLikelihood(contextT, likT);
  std::cout << "Initial likelihood is: " << phyloT->getValue() << std::endl;// calculate now because it will be needed for stochastic mapping
  
  auto params = phyloT->getSubstitutionModelParameters();
  std::cout <<"Parameters are:" << std::endl;
  for (size_t i = 0; i < params.size(); i++){
    std::cout <<"\t" << params[i].getName() << ": " << params[i].getValue() << std::endl;

  }
  std::vector<string> parameterNames = params.getParameterNames();
  vector<string> fixed_params_names;
  optimizeTraitModel(phyloT, 0.1, 5, parameterNames, fixed_params_names);
  std::cout << "final likelihood is: " << phyloT->getValue() << std::endl;
  delete traitAlpha;
  delete phyloT;
}
/******************************************************/
void verifyLikWeighted(){


  std::map<string, double> traitModelParams;
  traitModelParams["global_rate"] = 6.6919;
  std::string modelName = "singleRate";

  const IntegerAlphabet* traitAlpha = new IntegerAlphabet(1);
  vector<double> freqVals;
  double freq_0 = 0.438134;
  freqVals.push_back(freq_0);
  freqVals.push_back(1-freq_0);
  
  shared_ptr<IntegerFrequencySet> freqs = make_shared<FullIntegerFrequencySet>(traitAlpha, freqVals);
  std::shared_ptr<CharacterSubstitutionModel> model = setTraitModel(traitAlpha, freqs, modelName);
  setParametersToNewTraitModel(traitModelParams, model, freqs);
  auto& Q = model->getGenerator();
  MatrixTools::printForR(Q,"Q",std::cout);


  std::map<double, RowMatrix<double>> map_Pijt;
  vector<double> branchLengths = {0.1, 0.2, 0.3, 0.4};
  for (auto & branch: branchLengths){
    map_Pijt[branch] = model->getPij_t(branch);
  }
  double L45_0 = map_Pijt[0.3](0, 1)* map_Pijt[0.3](0, 0); // P01*Ls5 * P00Ls4
  double L45_1 = map_Pijt[0.3](1, 1)* map_Pijt[0.3](1, 0); // P11*Ls5 * P10Ls4
  double L12_0 = map_Pijt[0.1](0, 1)* map_Pijt[0.1](0, 1); // P01*Ls1 * P01Ls2
  double L12_1 = map_Pijt[0.1](1, 1)* map_Pijt[0.1](1, 1); // P11*Ls1 * P11Ls2

  double L123_0 = map_Pijt[0.4](0, 0) * ((map_Pijt[0.3](0,0) * L12_0) + (map_Pijt[0.3](0,1) * L12_1)); // (p00(0.3)L12(0) + p01(0.3)L12(1)) * p00(0.4)Ls3(0)
  double L123_1 = map_Pijt[0.4](1, 0) * ((map_Pijt[0.3](1,0) * L12_0) + (map_Pijt[0.3](1,1) * L12_1)); // (p10(0.3)L12(0) + p11(0.3)L12(1)) * p10(0.4)Ls3(0)

  // (p01(0.2)*L123_1 + p00(0.2)*L123_0) * (p00(0.3)*L45_0 + p01(0.3) * L45_1)
  double L12345_0 = ((map_Pijt[0.2](0, 1)*L123_1) + (map_Pijt[0.2](0, 0)*L123_0)) * ((map_Pijt[0.3](0, 1) * L45_1) + (map_Pijt[0.3](0, 0) * L45_0));
  // (p11(0.2)*L123_1 + p10(0.2)*L123_0) * (p11(0.3)*L45_1 + p10(0.3) * L45_0)
  double L12345_1 = ((map_Pijt[0.2](1, 1)*L123_1) + (map_Pijt[0.2](1, 0)*L123_0)) * ((map_Pijt[0.3](1, 1) * L45_1) + (map_Pijt[0.3](1, 0) * L45_0));

  double f_0 = L12345_0/(L12345_0 + L12345_1);
  double f_1 = 1 - f_0;


  std::cout << "F[0] = " << f_0 << std::endl;
  std::cout << "F[1] = " << f_1 << std::endl;
  //double f_0 = 0.3;
  //double f_1 = 0.7;
  //double f_0 = 1;
  //double f_1 = 0;

  auto& Q1 = model->getGenerator();
  MatrixTools::printForR(Q1,"Q",std::cout);

  double total_lik = (f_0 * L12345_0) + (f_1 * L12345_1);
  std::cout << "likelihood should be: " << total_lik << std::endl;
  std::cout << "log likelihood should be: " << std::log(total_lik) << std::endl;
  delete traitAlpha;

}
/******************************************************/

void verifyLikFixedFreqs(){

  std::map<string, double> traitModelParams;
  traitModelParams["global_rate"] = 1.41393;
  std::string modelName = "singleRate";

  const IntegerAlphabet* traitAlpha = new IntegerAlphabet(1);
  vector<double> freqVals;
  double freq_0 = 1e-08;
  freqVals.push_back(freq_0);
  freqVals.push_back(1-freq_0);
  
  shared_ptr<IntegerFrequencySet> freqs = make_shared<FullIntegerFrequencySet>(traitAlpha, freqVals);
  std::shared_ptr<CharacterSubstitutionModel> model = setTraitModel(traitAlpha, freqs, modelName);
  setParametersToNewTraitModel(traitModelParams, model, freqs);
  auto& Q = model->getGenerator();
  MatrixTools::printForR(Q,"Q",std::cout);


  std::map<double, RowMatrix<double>> map_Pijt;
  vector<double> branchLengths = {0.1, 0.2, 0.3, 0.4};
  for (auto & branch: branchLengths){
    map_Pijt[branch] = model->getPij_t(branch);
  }
  double L45_0 = map_Pijt[0.3](0, 1)* map_Pijt[0.3](0, 0); // P01*Ls5 * P00Ls4
  double L45_1 = map_Pijt[0.3](1, 1)* map_Pijt[0.3](1, 0); // P11*Ls5 * P10Ls4
  double L12_0 = map_Pijt[0.1](0, 1)* map_Pijt[0.1](0, 1); // P01*Ls1 * P01Ls2
  double L12_1 = map_Pijt[0.1](1, 1)* map_Pijt[0.1](1, 1); // P11*Ls1 * P11Ls2

  double L123_0 = map_Pijt[0.4](0, 0) * ((map_Pijt[0.3](0,0) * L12_0) + (map_Pijt[0.3](0,1) * L12_1)); // (p00(0.3)L12(0) + p01(0.3)L12(1)) * p00(0.4)Ls3(0)
  double L123_1 = map_Pijt[0.4](1, 0) * ((map_Pijt[0.3](1,0) * L12_0) + (map_Pijt[0.3](1,1) * L12_1)); // (p10(0.3)L12(0) + p11(0.3)L12(1)) * p10(0.4)Ls3(0)

  // (p01(0.2)*L123_1 + p00(0.2)*L123_0) * (p00(0.3)*L45_0 + p01(0.3) * L45_1)
  double L12345_0 = ((map_Pijt[0.2](0, 1)*L123_1) + (map_Pijt[0.2](0, 0)*L123_0)) * ((map_Pijt[0.3](0, 1) * L45_1) + (map_Pijt[0.3](0, 0) * L45_0));
  // (p11(0.2)*L123_1 + p10(0.2)*L123_0) * (p11(0.3)*L45_1 + p10(0.3) * L45_0)
  double L12345_1 = ((map_Pijt[0.2](1, 1)*L123_1) + (map_Pijt[0.2](1, 0)*L123_0)) * ((map_Pijt[0.3](1, 1) * L45_1) + (map_Pijt[0.3](1, 0) * L45_0));

  double f_0 = 1;
  double f_1 = 0;


  std::cout << "F[0] = " << f_0 << std::endl;
  std::cout << "F[1] = " << f_1 << std::endl;
  //double f_0 = 0.3;
  //double f_1 = 0.7;
  //double f_0 = 1;
  //double f_1 = 0;

  auto& Q1 = model->getGenerator();
  MatrixTools::printForR(Q1,"Q",std::cout);

  double total_lik = (f_0 * L12345_0) + (f_1 * L12345_1);
  std::cout << "likelihood should be: " << total_lik << std::endl;
  std::cout << "log likelihood should be: " << std::log(total_lik) << std::endl;
  std::cout << "fractional liks" << std::endl;
  auto lik_0 = f_0 * L12345_0;
  auto lik_1 = f_1 * L12345_1;
  std::cout << "log lik0 : " << std::log(lik_0) << std::endl;
  std::cout << "log lik1: " << std::log(lik_1) << std::endl;

  delete traitAlpha;

}
/*******************************************************/
double calculateSigma(vector<double> &u, vector<double> &V){
    double sigma = 0;
    size_t n = u.size();
    for (size_t i = 0; i < n; i++){
      double u_square = (u[i]*u[i]);
      //double V_square = (V[i]*V[i]);
      sigma += (u_square/V[i]);
    }
    return sigma/static_cast<double>(n);
}
double getExpectedTraitValue(double &t1, double &t2, double &y1, double &y2){
    double nominator = ((1/t1)*y1) + ((1/t2)*y2);
    double denominator = (1/t1) + (1/t2);
    return nominator/denominator;

}
double calculateNewBranchLength(double& originalBranchLength, double &t1, double &t2){
    double newBranchLength = originalBranchLength;
    newBranchLength += ((t1*t2)/(t1+t2));
    return newBranchLength;
}

void createNodeIdsTraitValuesMap(shared_ptr<PhyloTree> tree, std::unordered_map<uint, double> &traitValues, std::unordered_map<string, double> &leafTraitStates){
    auto leaves = tree->getAllLeaves();
    string leafName;
    for (auto &leaf : leaves){
        leafName = leaf->getName();
        traitValues[tree->getNodeIndex(leaf)] = leafTraitStates[leafName];
    }
}
void createBranchLengthsMap(shared_ptr<PhyloTree> tree, std::unordered_map<uint, double> &branchLengths){
    //auto father = tree_->getFatherOfNode(tree_->getNode(nodeId));
    //uint fatherId = tree_->getNodeIndex(father);
    auto nodeIds = tree->getNodeIndexes(tree->getAllNodes());
    for (auto &nodeId : nodeIds){
        if (nodeId == tree->getRootIndex()){
            continue;
        }
        if (tree->isLeaf(nodeId)){
            auto branch = tree->getEdgeToFather(nodeId);
            double branchLength = branch->getLength();
            branchLengths[nodeId] = branchLength;

        }else{
            branchLengths[nodeId] = -1; // should be updated during the recursion
        }

    }
}

void updateRequiredVectorsForLikCalculationRec(shared_ptr<PhyloTree> tree, uint nodeId, std::unordered_map<uint, double> &branchLengths, std::unordered_map<uint, double> &traitValues, vector<double> &u, vector<double> &V){
    if (tree->isLeaf(nodeId)){
        return;
    }
    auto sons = tree->getSons(nodeId);
    vector<double> sonsBranchLengths;
    for (auto &sonId : sons){
        if (!(tree->isLeaf(sonId))){
            updateRequiredVectorsForLikCalculationRec(tree, sonId, branchLengths, traitValues, u, V);
           
        }
    }
    double t1 = branchLengths[sons[0]];
    double t2 = branchLengths[sons[1]];
    if ((t1 < 0) || (t2 < 0)){
        throw Exception("updateRequiredVectorsForLikCalculationRec(): Something went wrong in the recursion!");

    }
    traitValues[nodeId] = getExpectedTraitValue(t1, t2, traitValues[sons[0]], traitValues[sons[1]]);
    u.push_back(traitValues[sons[0]]-traitValues[sons[1]]);
    V.push_back(t1+t2);
    if (nodeId == tree->getRootIndex()){
      double rootError = 0;
      branchLengths[nodeId] = calculateNewBranchLength(rootError, t1, t2);
    }else{    
      auto branch = tree->getEdgeToFather(nodeId);
      double branchLength = branch->getLength();
      branchLengths[nodeId] = calculateNewBranchLength(branchLength, t1, t2);

    }

}

double calculateLikelihood(shared_ptr<PhyloTree> tree, std::unordered_map<string, double> &leafTraitStates){
    std::unordered_map<uint, double> branchLengths;
    std::unordered_map<uint, double> traitValues;
    createNodeIdsTraitValuesMap(tree, traitValues, leafTraitStates);
    createBranchLengthsMap(tree, branchLengths);
    vector<double> V; // variances (same notation as in the paper)
    vector<double> u; // contrasts (same notation as in the paper)
    uint rootId = tree->getRootIndex();
    updateRequiredVectorsForLikCalculationRec(tree, rootId, branchLengths, traitValues, u, V);
    auto rootSons = tree->getSons(rootId);
    double rootError = 0;
    double rootUncertainty = calculateNewBranchLength(rootError, branchLengths[rootSons[0]], branchLengths[rootSons[1]]);
    u.push_back(0);
    V.push_back(rootUncertainty);
    double sigma = calculateSigma(u, V);
    std::cout << "sigma is " << sigma << std::endl;
    double likelihood = 0;
    for (size_t i = 0; i < u.size(); i++){
        likelihood += (std::log(V[i]*sigma) + ((u[i]*u[i])/(sigma * V[i])));
    }
    likelihood += static_cast<double>(u.size())*std::log(2* M_PI);
    likelihood *= (-0.5);
    // print the expected ancestral states:
    std::cout << "Expected ancestral states" << std::endl;
    auto it = traitValues.begin();
    while (it != traitValues.end()){
      std::cout << "node id " << it->first << " trait value: " << traitValues[it->first] << std::endl;
      it ++;
    }
    return likelihood;

}
/******************************************************/
void testBMTraitModel(){
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("((((P:0.21,H:0.21):0.28,M:0.49):0.13,A:0.62):0.38,G:1.0);", false, "", false, false));
  std::unordered_map<string, double> traitData;
  traitData["G"] = -1.46968;
  traitData["A"] = 2.02815;
  traitData["M"] = 2.37024;
  traitData["P"] = 3.61092;
  traitData["H"] = 4.09434;
  double likelihood = calculateLikelihood(pTree, traitData);
  std::cout << "likelihood is " << likelihood << std::endl;

}
/******************************************************/
void testBrownianMotionMLAncestralstates(){
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("((((P:0.21,H:0.21):0.28,M:0.49):0.13,A:0.62):0.38,G:1.0);", false, "", false, false));
  std::unordered_map<string, double> traitData;
  traitData["G"] = -1.46968;
  traitData["A"] = 2.02815;
  traitData["M"] = 2.37024;
  traitData["P"] = 3.61092;
  traitData["H"] = 4.09434;
  double muMLE = 1.18372461;
  double sigmaMLE = 3.153838667703223;
  // test the case where the parameters match the MLE
  BrownianMotionAncestralReconstruction* bmAncestral = new BrownianMotionAncestralReconstruction(muMLE, sigmaMLE, muMLE, sigmaMLE, traitData, pTree);
  bmAncestral->reconstructAncestralStates();
  auto reconstructedStates = bmAncestral->getAncestralStates();
  auto it = reconstructedStates.begin();
  while ( it != reconstructedStates.end()){
    std::cout << "Node id is " << it->first << " state is " << it->second << std::endl;
    it ++;
  }
  delete bmAncestral;

}
/******************************************************/
void testChromosomeBMSubstitutionModel(){
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((S1:0.1,S2:0.1):0.3,S3:0.4):0.2,(S4:0.3,S5:0.3):0.3);", false, "", false, false));
  const ChromosomeAlphabet* alphabetChr = new ChromosomeAlphabet(3,20);
  auto rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
  auto parTree =  std::make_shared<ParametrizablePhyloTree>(*pTree);
  std::map<int, std::vector<double>> mapOfParamValues;
  mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::GAIN)].push_back(2);
  mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::LOSS)].push_back(2);
  mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DUPL)].push_back(3);
  mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)].push_back(1.3);
  mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)].push_back(0.1);
  uint chrRange = 15;
  int baseNumber = 6;
  std::vector<int> rateChangeType;
  for (size_t i = 0; i < 5; i++){
    rateChangeType.push_back(static_cast<int>(ChromosomeNumberDependencyFunction::CONSTANT));

  }
  std::shared_ptr<VectorSiteContainer> vscChr = std::make_shared<VectorSiteContainer>(alphabetChr);
	vscChr->addSequence(BasicSequence("S1", "3", alphabetChr));
	vscChr->addSequence(BasicSequence("S2", "6", alphabetChr));
	vscChr->addSequence(BasicSequence("S3", "18", alphabetChr));
	vscChr->addSequence(BasicSequence("S4", "12", alphabetChr));
	vscChr->addSequence(BasicSequence("S5", "6", alphabetChr));
  std::shared_ptr<ChromosomeSubstitutionModel> chrModel1 = std::make_shared<ChromosomeSubstitutionModel>(alphabetChr, mapOfParamValues, baseNumber, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType, false);

  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
  subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdist->clone()), parTree);
  auto nodes = pTree->getAllNodes();
  std::vector<uint> modelNodes;
  for (size_t i = 0; i < nodes.size(); i++){
    auto nodeId = pTree->getNodeIndex(nodes[i]);
    if (nodeId == pTree->getRootIndex()){
      continue;
    }
    modelNodes.push_back(nodeId);
  }
  subProSim->addModel(chrModel1, modelNodes);
  SubstitutionProcess* subProcess= subProSim->clone();
  Context context;
  auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, *(vscChr->clone()), *subProcess, true);
  auto phylo1 = std::make_shared<SingleProcessPhyloLikelihood>(context, lik);
  double chrModelLik = phylo1->getValue();
  std::cout << "Regular chromosome model likelihood is: " << chrModelLik << std::endl;
  // defining chromosome model BM:
  double sigma = 3.15;
  double mu = 1.82;
  double traitState = 2;
  double minTraitState = 0.1;
  double maxTraitState = 10;
  std::map<int, std::vector<double>> mapOfParamValues2;
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::GAIN)].push_back(1);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::GAIN)].push_back(0.5);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::LOSS)].push_back(1.5);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::LOSS)].push_back(0.25);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::DUPL)].push_back(1.5);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::DUPL)].push_back(0.75);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)].push_back(1.3);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)].push_back(0.1);

  std::vector<int> rateChangeType2;
  rateChangeType2.push_back(static_cast<int>(ChromosomeNumberDependencyFunction::CONSTANT));
  for (size_t i = 0; i < 3; i++){
    rateChangeType2.push_back(static_cast<int>(ChromosomeNumberDependencyFunction::LINEAR));

  }
  rateChangeType2.push_back(static_cast<int>(ChromosomeNumberDependencyFunction::CONSTANT));
  std::shared_ptr<ChromosomeBMSubstitutionModel> chrBMModel = std::make_shared<ChromosomeBMSubstitutionModel>(mu, sigma, traitState, minTraitState, maxTraitState, minTraitState, maxTraitState, alphabetChr, mapOfParamValues2, baseNumber, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType2, false);
  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim2;
  auto parTree2 =  std::make_shared<ParametrizablePhyloTree>(*pTree);
  subProSim2 = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdist->clone()), parTree2);
  subProSim2->addModel(chrBMModel, modelNodes);
  SubstitutionProcess* subProcess2= subProSim2->clone();
  Context context2;
  auto lik2 = std::make_shared<LikelihoodCalculationSingleProcess>(context2, *(vscChr->clone()), *subProcess2, true);
  auto phylo2 = std::make_shared<SingleProcessPhyloLikelihood>(context2, lik2);
  double chrModelLik2 = phylo2->getValue();
  std::cout << "BM chromosome model likelihood is: " << chrModelLik2 << std::endl;

  //  printing the updated parameters

  auto parameters = phylo2->getSubstitutionModelParameters();
  for (size_t i = 0; i < parameters.size(); i++){
    std::cout << "\t" << parameters[i].getName() << " value is " << parameters[i].getValue() << std::endl;
  }

  phylo2->setParameterValue("Chromosome.state_1", 4);
  std::cout << "BM chromosome model likelihood is after changing state: " << phylo2->getValue() << std::endl;
  phylo1->setParameterValue("Chromosome.gain0_1", 3);
  phylo1->setParameterValue("Chromosome.loss0_1", 2.5);
  phylo1->setParameterValue("Chromosome.dupl0_1", 4.5);
  std::cout << "Regular chromosome model likelihood is after the parallel change: " << phylo1->getValue() << std::endl;
  // Now getting back to previous likelihood by changing the transition type parameters
  phylo2->setParameterValue("Chromosome.gain0_1", 1.5);
  phylo2->setParameterValue("Chromosome.gain1_1", 0.125);

  phylo2->setParameterValue("Chromosome.loss0_1", 1.8);
  phylo2->setParameterValue("Chromosome.loss1_1", 0.05);

  phylo2->setParameterValue("Chromosome.dupl0_1", 1);
  phylo2->setParameterValue("Chromosome.dupl1_1", 0.5);

  std::cout << "BM chromosome model likelihood is after changing back: " << phylo2->getValue() << std::endl;




}
void testBrownianMotionLikelihood(){
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("((((P:0.21,H:0.21):0.28,M:0.49):0.13,A:0.62):0.38,G:1.0);", false, "", false, false));
  std::unordered_map<string, double> traitData;
  traitData["G"] = -1.46968;
  traitData["A"] = 2.02815;
  traitData["M"] = 2.37024;
  traitData["P"] = 3.61092;
  traitData["H"] = 4.09434;
  double muMLE = 1.18372461;
  double sigmaMLE = 3.153838667703223;
  BrownianMotionLikelihood* bm = new BrownianMotionLikelihood(pTree, traitData);
  //bm->calculateContrasts();
  std::cout << "mu is: " << bm->getMuMLE() << std::endl;
  std::cout << "sigma is: " << bm->getSigmaMLE() << std::endl;
  bm->calculateLikelihood(muMLE, sigmaMLE);
  std::cout << "likelihood is: " << bm->getLikelihood() << std::endl;
  delete bm;

}
/******************************************************/
int main(){
  //testTraitWeighted();
  //testBrownianMotionMLAncestralstates();
  //testBrownianMotionLikelihood();
  testChromosomeBMSubstitutionModel();

  return 0;
  // string simulations_path = "/home/anat/Docs/Sida/simulations";
  // size_t num_of_simulations = 100;
  // double totalNumOfTransitions = 0;
  // std::map<string, double> state_changes;
  // state_changes["gain"] = 0;
  // state_changes["loss"] = 0;
  // state_changes["dupl"] = 0;
  // std::regex states_pattern("from state: (\\d+)[^\\d]+t = [\\d.]+[^\\d]+to state = (\\d+)");
  // for (size_t i = 0; i < num_of_simulations; i++){
  //   string evolPath = simulations_path + "/" + std::to_string(i) + '/' + "simulatedEvolutionPaths.txt";
  //   std::ifstream file(evolPath); // Open the file

  //   if (!file.is_open()) {
  //       std::cerr << "Error opening file: " << evolPath << std::endl;
  //       return 1;
  //   }

  //   // Read the entire file into a string
  //   std::string content((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));
  //   file.close(); 
  //   std::regex pattern("Total number of transitions is: (\\d+)");
  //   std::smatch matches;
  //   if (std::regex_search(content, matches, pattern)) {
  //       // Print the matched sentence
  //       int numOfTransitions = std::stoi(matches[1]);
  //       totalNumOfTransitions += numOfTransitions;
  //   } else {
  //       std::cerr << "Sentence not found in the text." << std::endl;
  //       return 1;
  //   }

  //   // Search for matches in the text
  //   std::smatch matches_states_test;
  //   if (!(std::regex_search(content, matches_states_test, states_pattern))) {
  //     std::cerr << "Sentence not found in the text." << std::endl;
  //   }
  //   std::smatch matches_states;
  //   std::string::const_iterator searchStart(content.cbegin());
  //   while (std::regex_search(searchStart, content.cend(), matches_states, states_pattern)) {
  //       // Extract the numbers
  //       int fromState = std::stoi(matches_states[1]);
  //       int toState = std::stoi(matches_states[2]);
  //       if (toState == fromState + 1){
  //         state_changes["gain"] += 1;
  //       }else if (toState == fromState * 2){
  //         state_changes["dupl"] += 1;
  //       }else if (toState == fromState-1){
  //         state_changes["loss"] += 1;
  //       }else{
  //         searchStart = matches_states.suffix().first;
  //         continue;
          
  //       }
  //       // Update the searchStart iterator to search for the next match
  //       searchStart = matches_states.suffix().first;
  //   }
  // }
  // totalNumOfTransitions /= num_of_simulations;
  // std::cout << "Average number of transitions " << totalNumOfTransitions << std::endl;
  // auto it = state_changes.begin();
  // while ( it != state_changes.end()){
  //   state_changes[it->first] /= (double)num_of_simulations;
  //   std::cout << it->first << " : " << it->second << std::endl;
  //   it ++;
  // }

  // return 0;
}
  //fix seed for debugging purposes
  // double seedUb = 10000000;
  // RandomTools::setSeed(static_cast<long int>(seedUb));
  // // defining the models
  // Newick reader;
  // shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((S1:0.1,S2:0.1):0.3,S3:0.4):0.2,(S4:0.3,S5:0.3):0.3);", false, "", false, false));
  // std::map<string, double> traitModelParams;
  // traitModelParams["global_rate"] = 2.0;
  // testMulti1(traitModelParams, pTree);
  
  // test chi square
  // double chi = RandomTools::qChisq(0.95, 20);;
  // std::cout << "chi square" << std::endl;
  // std::cout << chi << std::endl;
  // const BinaryAlphabet* alphabet = new BinaryAlphabet();
  // double mu = 1.;
  // double pi0 = 0.25;
  // auto twoParamModel = std::make_shared<TwoParameterBinarySubstitutionModel>(alphabet,mu,pi0);
  // std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
  // vector <double> rootFreqs;
  // rootFreqs.push_back(0.25);
  // rootFreqs.push_back(0.75);
  // std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(twoParamModel->getStateMap(), false)), rootFreqs);
  // std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
  // std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdist->clone()), parTree, rootFrequencies);        

  // // process character data
  // std::shared_ptr<VectorSiteContainer> sites = std::make_shared<VectorSiteContainer>(alphabet);
	// sites->addSequence(BasicSequence("S1", "1", alphabet));
	// sites->addSequence(BasicSequence("S2", "1", alphabet));
	// sites->addSequence(BasicSequence("S3", "0", alphabet));
	// sites->addSequence(BasicSequence("S4", "0", alphabet));
	// sites->addSequence(BasicSequence("S5", "1", alphabet));

  // auto nodes = pTree->getAllNodes();
  // std::vector<uint> modelNodes;
  // for (size_t i = 0; i < nodes.size(); i++){
  //   auto nodeId = pTree->getNodeIndex(nodes[i]);
  //   if (nodeId == pTree->getRootIndex()){
  //     continue;
  //   }
  //   modelNodes.push_back(nodeId);
  // }
  // subProSim->addModel(std::shared_ptr<TwoParameterBinarySubstitutionModel>(twoParamModel->clone()), modelNodes);
		        
  // // create tree likelihood function for the trait model
  // Context context1;
  // auto nsubPro1 = subProSim->clone();
  // auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context1, *sites->clone(), *nsubPro1);
  // auto phylo1 = std::make_shared<SingleProcessPhyloLikelihood>(context1, lik);

  // phylo1->getValue();

  // // get expected tree for the chromosome model
  // StochasticMapping* stm = new StochasticMapping(lik, 3);
  // auto treeChr = createExpectedMappingHistory(lik, 3, stm);
  // // creating the chromosome number model
  // const ChromosomeAlphabet* alphabetChr = new ChromosomeAlphabet(3,20);
  // auto rdist2 = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
  // auto parTree2 =  std::make_shared<ParametrizablePhyloTree>(*treeChr);
  // std::map<int, std::vector<double>> mapOfParamValues;
  // mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::GAIN)].push_back(2);
  // mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::LOSS)].push_back(2);
  // mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DUPL)].push_back(3);
  // mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)].push_back(1.3);
  // mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)].push_back(0.1);
  // uint chrRange = 15;
  // int baseNumber = 6;
  // std::vector<int> rateChangeType;
  // for (size_t i = 0; i < 5; i++){
  //   rateChangeType.push_back(static_cast<int>(ChromosomeNumberDependencyFunction::CONSTANT));

  // }

  // std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim2;
  // std::shared_ptr<ChromosomeSubstitutionModel> chrModel1 = std::make_shared<ChromosomeSubstitutionModel>(alphabetChr, mapOfParamValues, baseNumber, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType);
  // std::map<int, std::vector<double>> mapOfParamValues2;
  // mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::GAIN)].push_back(1);
  // mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::LOSS)].push_back(1);
  // mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::DUPL)].push_back(0.1);
  // mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)].push_back(0.1);
  // mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)].push_back(0.1);
  // std::shared_ptr<ChromosomeSubstitutionModel> chrModel2 = std::make_shared<ChromosomeSubstitutionModel>(alphabetChr, mapOfParamValues2, baseNumber, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType);

  // subProSim2 = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdist2->clone()), parTree2);
  // auto nodesMap = getNodesForEachModel(treeChr, stm);
  // subProSim2->addModel(chrModel1, nodesMap[0]);
  // subProSim2->addModel(chrModel2, nodesMap[1]);

  // SubstitutionProcess* nsubPro2= subProSim2->clone();
  // Context context2;
  // // chromosome number data
  // std::shared_ptr<VectorSiteContainer> sites2 = std::make_shared<VectorSiteContainer>(alphabetChr);
	// sites2->addSequence(BasicSequence("S1-1", "3", alphabetChr));
	// sites2->addSequence(BasicSequence("S2-1", "6", alphabetChr));
	// sites2->addSequence(BasicSequence("S3-0", "18", alphabetChr));
	// sites2->addSequence(BasicSequence("S4-0", "12", alphabetChr));
	// sites2->addSequence(BasicSequence("S5-1", "6", alphabetChr));

  // auto lik2 = std::make_shared<LikelihoodCalculationSingleProcess>(context2, *sites2->clone(), *nsubPro2, true);
  // auto phylo2 = std::make_shared<SingleProcessPhyloLikelihood>(context2, lik2, lik2->getParameters());
  // std::cout << "Likelihoods separately:" << std::endl;

  // auto lik_val1 = phylo1->getValue();
  // auto lik_val2 = phylo2->getValue();

  // std::cout << "Trait likelihood: " << lik_val1 << std::endl;
  // std::cout << "Chromosome number likelihood: " << lik_val2 << std::endl;

  // // Creating model collection instance
  // //std::shared_ptr<const FrequencySet> rootFreqs2 = lik2->getSubstitutionProcess().getRootFrequencySet();

  // auto modelColl=std::make_shared<SubstitutionProcessCollection>();
  
  // modelColl->addModel(twoParamModel, 1);
  // modelColl->addModel(chrModel1, 2);
  // modelColl->addModel(chrModel2, 3);

  // modelColl->addFrequencies(rootFrequencies, 1);
  // //auto root_freq_cloned = rootFreqs2->clone();
  // //auto rootFreqs_chr = std::shared_ptr<FrequencySet>(root_freq_cloned);
  // //modelColl->addFrequencies(rootFreqs_chr, 2);
  // modelColl->addDistribution(rdist, 1);
  // modelColl->addDistribution(rdist2, 2);

  // modelColl->addTree(parTree, 1);
  // modelColl->addTree(parTree2, 2);
  // vector<uint> vP1m1 = modelNodes;

  // map<size_t, Vuint> mModBr1;
  // mModBr1[1]=vP1m1;
 

  // modelColl->addSubstitutionProcess(1, mModBr1, 1, 1);
                                   
  // map<size_t, Vuint> mModBr2;
  // std::vector<uint> vP2m1 = nodesMap[0];
  // std::vector<uint> vP2m2 = nodesMap[1];
  // mModBr2[2]=vP2m1;
  // mModBr2[3]=vP2m2;

  // modelColl->addSubstitutionProcess(2, mModBr2, 2, 2);

  // // Now setting the joint likelihood object
  //   // Likelihoods
  // Context context;
  // SubstitutionProcess* sP1c=nsubPro1->clone();
  // SubstitutionProcess* sP2c=nsubPro2->clone();
  // auto pc(std::make_shared<PhyloLikelihoodContainer>(context, *modelColl));



  // auto lik1_j = std::make_shared<LikelihoodCalculationSingleProcess>(context, *(sites->clone()), *sP1c);

  // pc->addPhyloLikelihood(1, new SingleProcessPhyloLikelihood(context, lik1_j));
    
  // auto lik2_j = std::make_shared<LikelihoodCalculationSingleProcess>(context, *(sites2->clone()), *sP2c, true);

  // pc->addPhyloLikelihood(2, new SingleProcessPhyloLikelihood(context, lik2_j));
  // //JointTraitChromosomeLikelihood tl(context, pc, true, true, 1000, "Binary", rateChangeType, true, "Ranges", sites2->clone());

  // //auto jl = tl.getValue();
  // //std::cout << "Joint likelihood is " << jl << std::endl;
  // std::cout << "sum of likelihoods is: " << lik_val1 + lik_val2 << std::endl;
  // //tl.optimize(0.01, 5, 5);


  // // tl.setParamUpdateMode(false);
  // // tl.setParameterValue("Chromosome.baseNum_1", 3);
  // // tl.setParameterValue("Chromosome.dupl0_1", 0.01);
  // // std::cout << "After parameter values update in chromosome model " << tl.getValue() << std::endl;
  // // phylo2->setParameterValue("Chromosome.baseNum_1", 3);
  // // phylo2->setParameterValue("Chromosome.dupl0_1", 0.01);
  // // std::cout << "chromosome number update separate chromosome model" << phylo2->getValue() << std::endl;
  // // std::cout << "chromosome number update separate trait model " << phylo1->getValue() << std::endl;
  // // std::cout << "sum of likelihoods is " << phylo2->getValue() + phylo1->getValue() << std::endl;

  // // tl.setParamUpdateMode(true);
  // // tl.setParameterValue("TwoParameterBinary.mu_1_1", 0.25);
  // // tl.setParameterValue("TwoParameterBinary.pi0_1_1", 0.8);

  // // std::cout << "After parameter values update " << tl.getValue() << std::endl;
  // // auto lik_val_T =  tl.getAbstractPhyloLikelihood(tl.getNumbersOfPhyloLikelihoods()[0])->getValue();
  // // auto lik_val_C = tl.getPhylo2()->getValue();
  // // std::cout << "trait lik !!! " << lik_val_T << std::endl;
  // // std::cout << "chromosome lik !!! " << lik_val_C << std::endl;
  // // std::cout << "overall lik: "  << (lik_val_T + lik_val_C) << std::endl;
  // // tl.setParamUpdateMode(false);
  // // tl.setParameterValue("Chromosome.dupl0_1", 7);
  // // std::cout << "After parameter values update 2!" << tl.getValue() << std::endl;
  // // auto lik_val_T2 =  tl.getAbstractPhyloLikelihood(tl.getNumbersOfPhyloLikelihoods()[0])->getValue();
  // // auto lik_val_C2 = tl.getPhylo2()->getValue();
  // // std::cout << "trait lik 2!!! " << lik_val_T2 << std::endl;
  // // std::cout << "chromosome lik 2!!! " << lik_val_C2 << std::endl;
  // // std::cout << "overall lik2: "  << (lik_val_T2 + lik_val_C2) << std::endl;


  


  // // // get likelihood using separate phylo-likelihoods
  // // phylo1->setParameterValue("TwoParameterBinary.mu_1", 0.25);
  // // phylo1->setParameterValue("TwoParameterBinary.pi0_1", 0.8);
  

  // // auto lik_val_trait = phylo1->getValue();
  // // std::cout << "Likelihood only of the trait model: " << lik_val_trait << std::endl;
  
  // // // printing the updated parameters

  // // auto parameters = tl.getParameters();
  // // auto paramsNames = parameters.getParameterNames();
  // // for (size_t i = 0; i < paramsNames.size(); i++){
  // //   std::cout << paramsNames[i] << " " << tl.getParameter(paramsNames[i]).getValue() << std::endl;
  // // }


  // // // now to see whether the parameters are updated in each iteration
  // // auto likTrait = std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(tl.getAbstractPhyloLikelihood(tl.getNumbersOfPhyloLikelihoods()[0])->getLikelihoodCalculation());
  // // parameters = likTrait->getParameters();
  // // std::cout << "TwoParameterBinary.mu_1 " << likTrait->getParameter("TwoParameterBinary.mu_1").getValue() << std::endl;
  // // std::cout << "TwoParameterBinary.pi0_1 " << likTrait->getParameter("TwoParameterBinary.pi0_1").getValue() << std::endl;
  






  
  //JointPhyloLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool expectedHistory, bool weightedFrequencies, size_t numOfMappings,
