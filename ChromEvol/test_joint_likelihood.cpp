#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/JointPhyloLikelihood.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include "ChromosomeNumberOptimizer.h"
#include "JointTraitChromosomeLikelihood.h"

using namespace bpp;
using namespace std;


std::shared_ptr<PhyloTree> createExpectedMappingHistory(std::shared_ptr<LikelihoodCalculationSingleProcess> lik, uint mappingsNum, StochasticMapping* stm){
  
    //ChromEvolOptions::NumOfSimulations_);
    stm->generateStochasticMapping();
    std::vector<std::shared_ptr<PhyloTree>> mappings;
    for (size_t i = 0; i < mappingsNum; i++){
        auto mappingTree = stm->createMappingHistoryTree(i);
        mappings.push_back(mappingTree);
        std::cout << "****************" << std::endl;
    }
    std::cout << "Expected mapping tree: " << std::endl;
        
    auto expected_mapping = stm->generateExpectedMapping(mappings);
    return expected_mapping;
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
        nodeModels[static_cast<uint>(nodeState)].push_back(sons[j]);
      }
    }
  }
  return nodeModels;
}

int main(){
  //fix seed for debugging purposes
  double seedUb = 10000000;
  RandomTools::setSeed(static_cast<long int>(seedUb));
  // defining the models
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((S1:0.1,S2:0.1):0.3,S3:0.4):0.2,(S4:0.3,S5:0.3):0.3);", false, "", false, false));
  std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*pTree);
        
  // create a binary model
  const BinaryAlphabet* alphabet = new BinaryAlphabet();
  double mu = 1.;
  double pi0 = 0.5;
  auto twoParamModel = std::make_shared<TwoParameterBinarySubstitutionModel>(alphabet,mu,pi0);
  std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
  vector <double> rootFreqs;
  rootFreqs.push_back(0.25);
  rootFreqs.push_back(0.75);
  std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(twoParamModel->getStateMap(), false)), rootFreqs);
  std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdist->clone()), parTree, rootFrequencies);        

  // process character data
  std::shared_ptr<VectorSiteContainer> sites = std::make_shared<VectorSiteContainer>(alphabet);
	sites->addSequence(BasicSequence("S1", "1", alphabet));
	sites->addSequence(BasicSequence("S2", "1", alphabet));
	sites->addSequence(BasicSequence("S3", "0", alphabet));
	sites->addSequence(BasicSequence("S4", "0", alphabet));
	sites->addSequence(BasicSequence("S5", "1", alphabet));

  auto nodes = pTree->getAllNodes();
  std::vector<uint> modelNodes;
  for (size_t i = 0; i < nodes.size(); i++){
    auto nodeId = pTree->getNodeIndex(nodes[i]);
    if (nodeId == pTree->getRootIndex()){
      continue;
    }
    modelNodes.push_back(nodeId);
  }
  subProSim->addModel(std::shared_ptr<TwoParameterBinarySubstitutionModel>(twoParamModel->clone()), modelNodes);
		        
  // create tree likelihood function for the trait model
  Context context1;
  auto nsubPro1 = subProSim->clone();
  auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context1, *sites->clone(), *nsubPro1);
  auto phylo1 = std::make_shared<SingleProcessPhyloLikelihood>(context1, lik);

  phylo1->getValue();

  // get expected tree for the chromosome model
  StochasticMapping* stm = new StochasticMapping(lik, 3);
  auto treeChr = createExpectedMappingHistory(lik, 3, stm);
  // creating the chromosome number model
  const ChromosomeAlphabet* alphabetChr = new ChromosomeAlphabet(3,20);
  auto rdist2 = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
  auto parTree2 =  std::make_shared<ParametrizablePhyloTree>(*treeChr);
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

  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim2;
  std::shared_ptr<ChromosomeSubstitutionModel> chrModel1 = std::make_shared<ChromosomeSubstitutionModel>(alphabetChr, mapOfParamValues, baseNumber, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType);
  std::map<int, std::vector<double>> mapOfParamValues2;
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::GAIN)].push_back(1);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::LOSS)].push_back(1);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::DUPL)].push_back(0.1);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)].push_back(0.1);
  mapOfParamValues2[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)].push_back(0.1);
  std::shared_ptr<ChromosomeSubstitutionModel> chrModel2 = std::make_shared<ChromosomeSubstitutionModel>(alphabetChr, mapOfParamValues2, baseNumber, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType);

  subProSim2 = std::make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdist2->clone()), parTree2);
  auto nodesMap = getNodesForEachModel(treeChr, stm);
  subProSim2->addModel(chrModel1, nodesMap[0]);
  subProSim2->addModel(chrModel2, nodesMap[1]);

  SubstitutionProcess* nsubPro2= subProSim2->clone();
  Context context2;
  // chromosome number data
  std::shared_ptr<VectorSiteContainer> sites2 = std::make_shared<VectorSiteContainer>(alphabetChr);
	sites2->addSequence(BasicSequence("S1-1", "3", alphabetChr));
	sites2->addSequence(BasicSequence("S2-1", "6", alphabetChr));
	sites2->addSequence(BasicSequence("S3-0", "18", alphabetChr));
	sites2->addSequence(BasicSequence("S4-0", "12", alphabetChr));
	sites2->addSequence(BasicSequence("S5-1", "6", alphabetChr));

  auto lik2 = std::make_shared<LikelihoodCalculationSingleProcess>(context2, *sites2->clone(), *nsubPro2, true);
  auto phylo2 = std::make_shared<SingleProcessPhyloLikelihood>(context2, lik2, lik2->getParameters());
  std::cout << "Likelihoods separately:" << std::endl;

  auto lik_val1 = phylo1->getValue();
  auto lik_val2 = phylo2->getValue();

  std::cout << "Trait likelihood: " << lik_val1 << std::endl;
  std::cout << "Chromosome number likelihood: " << lik_val2 << std::endl;

  // Creating model collection instance
  //std::shared_ptr<const FrequencySet> rootFreqs2 = lik2->getSubstitutionProcess().getRootFrequencySet();

  auto modelColl=std::make_shared<SubstitutionProcessCollection>();
  
  modelColl->addModel(twoParamModel, 1);
  modelColl->addModel(chrModel1, 2);
  modelColl->addModel(chrModel2, 3);

  modelColl->addFrequencies(rootFrequencies, 1);
  //auto root_freq_cloned = rootFreqs2->clone();
  //auto rootFreqs_chr = std::shared_ptr<FrequencySet>(root_freq_cloned);
  //modelColl->addFrequencies(rootFreqs_chr, 2);
  modelColl->addDistribution(rdist, 1);
  modelColl->addDistribution(rdist2, 2);

  modelColl->addTree(parTree, 1);
  modelColl->addTree(parTree2, 2);
  vector<uint> vP1m1 = modelNodes;

  map<size_t, Vuint> mModBr1;
  mModBr1[1]=vP1m1;
 

  modelColl->addSubstitutionProcess(1, mModBr1, 1, 1);
                                   
  map<size_t, Vuint> mModBr2;
  std::vector<uint> vP2m1 = nodesMap[0];
  std::vector<uint> vP2m2 = nodesMap[1];
  mModBr2[2]=vP2m1;
  mModBr2[3]=vP2m2;

  modelColl->addSubstitutionProcess(2, mModBr2, 2, 2);

  // Now setting the joint likelihood object
    // Likelihoods
  Context context;
  SubstitutionProcess* sP1c=nsubPro1->clone();
  SubstitutionProcess* sP2c=nsubPro2->clone();
  auto pc(std::make_shared<PhyloLikelihoodContainer>(context, *modelColl));



  auto lik1_j = std::make_shared<LikelihoodCalculationSingleProcess>(context, *(sites->clone()), *sP1c);

  pc->addPhyloLikelihood(1, new SingleProcessPhyloLikelihood(context, lik1_j));
    
  auto lik2_j = std::make_shared<LikelihoodCalculationSingleProcess>(context, *(sites2->clone()), *sP2c, true);

  pc->addPhyloLikelihood(2, new SingleProcessPhyloLikelihood(context, lik2_j));
  //JointTraitChromosomeLikelihood tl(context, pc, true, true, 1000, "Binary", rateChangeType, true, "Ranges", sites2->clone());

  //auto jl = tl.getValue();
  //std::cout << "Joint likelihood is " << jl << std::endl;
  std::cout << "sum of likelihoods is: " << lik_val1 + lik_val2 << std::endl;
  //tl.optimize(0.01, 5, 5);


  // tl.setParamUpdateMode(false);
  // tl.setParameterValue("Chromosome.baseNum_1", 3);
  // tl.setParameterValue("Chromosome.dupl0_1", 0.01);
  // std::cout << "After parameter values update in chromosome model " << tl.getValue() << std::endl;
  // phylo2->setParameterValue("Chromosome.baseNum_1", 3);
  // phylo2->setParameterValue("Chromosome.dupl0_1", 0.01);
  // std::cout << "chromosome number update separate chromosome model" << phylo2->getValue() << std::endl;
  // std::cout << "chromosome number update separate trait model " << phylo1->getValue() << std::endl;
  // std::cout << "sum of likelihoods is " << phylo2->getValue() + phylo1->getValue() << std::endl;

  // tl.setParamUpdateMode(true);
  // tl.setParameterValue("TwoParameterBinary.mu_1_1", 0.25);
  // tl.setParameterValue("TwoParameterBinary.pi0_1_1", 0.8);

  // std::cout << "After parameter values update " << tl.getValue() << std::endl;
  // auto lik_val_T =  tl.getAbstractPhyloLikelihood(tl.getNumbersOfPhyloLikelihoods()[0])->getValue();
  // auto lik_val_C = tl.getPhylo2()->getValue();
  // std::cout << "trait lik !!! " << lik_val_T << std::endl;
  // std::cout << "chromosome lik !!! " << lik_val_C << std::endl;
  // std::cout << "overall lik: "  << (lik_val_T + lik_val_C) << std::endl;
  // tl.setParamUpdateMode(false);
  // tl.setParameterValue("Chromosome.dupl0_1", 7);
  // std::cout << "After parameter values update 2!" << tl.getValue() << std::endl;
  // auto lik_val_T2 =  tl.getAbstractPhyloLikelihood(tl.getNumbersOfPhyloLikelihoods()[0])->getValue();
  // auto lik_val_C2 = tl.getPhylo2()->getValue();
  // std::cout << "trait lik 2!!! " << lik_val_T2 << std::endl;
  // std::cout << "chromosome lik 2!!! " << lik_val_C2 << std::endl;
  // std::cout << "overall lik2: "  << (lik_val_T2 + lik_val_C2) << std::endl;


  


  // // get likelihood using separate phylo-likelihoods
  // phylo1->setParameterValue("TwoParameterBinary.mu_1", 0.25);
  // phylo1->setParameterValue("TwoParameterBinary.pi0_1", 0.8);
  

  // auto lik_val_trait = phylo1->getValue();
  // std::cout << "Likelihood only of the trait model: " << lik_val_trait << std::endl;
  
  // // printing the updated parameters

  // auto parameters = tl.getParameters();
  // auto paramsNames = parameters.getParameterNames();
  // for (size_t i = 0; i < paramsNames.size(); i++){
  //   std::cout << paramsNames[i] << " " << tl.getParameter(paramsNames[i]).getValue() << std::endl;
  // }


  // // now to see whether the parameters are updated in each iteration
  // auto likTrait = std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(tl.getAbstractPhyloLikelihood(tl.getNumbersOfPhyloLikelihoods()[0])->getLikelihoodCalculation());
  // parameters = likTrait->getParameters();
  // std::cout << "TwoParameterBinary.mu_1 " << likTrait->getParameter("TwoParameterBinary.mu_1").getValue() << std::endl;
  // std::cout << "TwoParameterBinary.pi0_1 " << likTrait->getParameter("TwoParameterBinary.pi0_1").getValue() << std::endl;
  






  
  //JointPhyloLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool expectedHistory, bool weightedFrequencies, size_t numOfMappings,
  return 0;
}
