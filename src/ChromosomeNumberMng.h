//
// File: ChromosomeNumberMng.h
// Created by: Anat Shafir
// Created on: Mon September 11 14:57 2020
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/
#ifndef CHROMEVOL_CHROMOSOMENUMBERMNG_H
#define CHROMEVOL_CHROMOSOMENUMBERMNG_H
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>



//from bpp-seq
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Io/AbstractISequence.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io/Pasta.h>

//from bpp-phyl
#include "LikelihoodUtils.h"
#include "StochasticMappingUtils.h"
#include <Bpp/Phyl/Tree/PhyloTreeTools.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/MarginalAncestralReconstruction.h>
#include <Bpp/Phyl/Likelihood/JointMLAncestralReconstruction.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowNumeric.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>
#include "ChromosomeNumberOptimizer.h"
#include "ComputeChromosomeTransitionsExp.h"
#include "ChromosomeSubstitutionModel.h"
#include "ChrFasta.h"
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSequenceSimulator.h>
#include "ChromosomeTraitOptimizer.h"
#include "JointChromosomeBMMng.h"
#include "BrownianMotionSimulator.h"



//standard libraries
#include <string>
#include <vector>
#include <iostream>
#include <time.h>
#include <sys/stat.h>
#include <regex>
#include <experimental/filesystem>
#include <type_traits>
#include <set>

using namespace std;

namespace bpp{
    class ChromosomeNumberMng{
        private:
            PhyloTree* tree_;
            ChromosomeAlphabet* alphabet_;
            VectorSiteContainer* vsc_;
            std::map<uint, uint> chrRange_; //maxObserved-minObserved chromosome number
            unsigned int numberOfUniqueStates_; // number of unique states (number of chromosomes)
            bool nullHypothesisRejected_;
        public :
            typedef Table<double> DataTable;




        public:
            //constructor
            ChromosomeNumberMng(): tree_(0), alphabet_(0), vsc_(0), chrRange_(), numberOfUniqueStates_(0), nullHypothesisRejected_(false){}
            ChromosomeNumberMng(const ChromosomeNumberMng& mng):
                tree_(mng.tree_->clone()), alphabet_(mng.alphabet_->clone()), vsc_(mng.vsc_->clone()), chrRange_(mng.chrRange_), numberOfUniqueStates_(mng.numberOfUniqueStates_)
            {}
            ChromosomeNumberMng& operator=(const ChromosomeNumberMng& mng){
                tree_ = tree_->clone();
                alphabet_ = mng.alphabet_->clone();
                vsc_ = mng.vsc_->clone();
                chrRange_ = mng.chrRange_;
                numberOfUniqueStates_ = mng.numberOfUniqueStates_;
                nullHypothesisRejected_ = mng.nullHypothesisRejected_;
                return *this;
            }
            virtual ~ChromosomeNumberMng(){
                if (vsc_){
                    delete vsc_;
                }
                delete tree_;
                delete alphabet_;
            };

            //Functions for initialization of the model
            void getCharacterData(const string &path);
            VectorSiteContainer* getTraitData(const string &path=ChromEvolOptions::traitFilePath_) const;
            static void setMaxChrNum(unsigned int maxNumberOfChr);
            static void setMinChrNum(unsigned int minNumberOfChr);
            void getTree(const string &path, double treeLength = 0);

            // getters for testers
            const ChromosomeAlphabet* getAlphabet() const {return alphabet_;}
            const VectorSiteContainer* getSeqData() const {return vsc_;}
            const map<uint, uint> getChromosomeRange() const {return chrRange_;}
            const PhyloTree* getPhyloTree() const {return tree_;}
            


            //core functions of ChromEvol
            //void runTest();
            void runChromEvol();
            void printStochasticMappingResults(StochasticMapping* stm, Vdouble &dwellingTimesPerState, std::map<pair<size_t, size_t>, double> &numOfOccurencesPerTransition, VVdouble &ratesPerTransition, std::map<int, double> &expectationsTotal, const string &outStMappingPath);
            void runStochasticMapping(ChromosomeNumberOptimizer* chrOptimizer);
            ChromosomeNumberOptimizer* optimizeLikelihoodMultiStartPoints(bool brownianMotionNull = false) const;
            void runJointTraitChromosomeAnalysis();
            JointChromosomeBMMng* runJointChromosomeBMModel();
            void runJointChromosomeContinuousBMAnalysis();
            void getJointMLAncestralReconstruction(ChromosomeNumberOptimizer* optimizer, int* inferredRootState, VectorSiteContainer* vsc, ChromosomeTraitOptimizer* traitOpt=0) const;
            void getMarginalAncestralReconstruction(ChromosomeNumberOptimizer* chrOptimizer, const string &filePath, ChromosomeTraitOptimizer* traitOpt, JointChromosomeBMMng* bmMng = 0);
            void computeExpectations(ChromosomeNumberOptimizer* chrOptimizer, int numOfSimulations, VectorSiteContainer* chromsomeVsc, ChromosomeTraitOptimizer* traitOpt, JointChromosomeBMMng* bmMng = 0);
            void simulateData(string &characterPath);
            //whichSimulation parameter is by deault 0, which means that a single type of data should be simulated. However, when a joint trait model is
            // simulated we need to types of data - 1 refers to trait data, and 2 refers to chromosome data.
            void simulateData(bool into_dirs, size_t simNum, size_t &count_failed, SimpleSubstitutionProcessSiteSimulator* simulator, Alphabet* alpha, SiteSimulationResult** simResult, size_t whichSimulation=0);
            void simulateTraitData(SimpleSubstitutionProcessSiteSimulator** simulator, Alphabet** alpha, std::shared_ptr<ParametrizablePhyloTree> parTree, std::shared_ptr<FrequencySet> &rootFrequencies, std::shared_ptr<NonHomogeneousSubstitutionProcess> &subProSim);
            bool simulateChromosomeData(SimpleSubstitutionProcessSiteSimulator** simulator, Alphabet** alpha, std::shared_ptr<ParametrizablePhyloTree> parTree, std::shared_ptr<FrequencySet> &rootFrequencies, std::shared_ptr<NonHomogeneousSubstitutionProcess> &subProSim, std::shared_ptr<BrownianMotionSimulator> bmSim = nullptr);

            void printSimulatedData(vector<size_t> leavesStates, vector<string> leavesNames, size_t iter, string &countsPath, Alphabet* alpha, string prefix);
            void printTreeWithStates(PhyloTree tree, std::map<uint, std::vector<size_t>> &ancestors, const string &filePath, string prefix, std::map<uint, string>* prevNames = 0) const;
            void convertNodesNames(PhyloTree &tree, uint nodeId, std::map<uint, std::vector<size_t>> &ancestors, string prefix, bool alphabetStates = true, std::map<uint, string>* prevNames = 0) const;
            void writeOutputToFile(ChromosomeNumberOptimizer* chrOptimizer, int &inferrredRootState) const;
            void printLikParameters(ChromosomeNumberOptimizer* chrOptimizer, SingleProcessPhyloLikelihood* lik, ofstream &outFile) const;
            uint findMinCladeSize(std::map<uint, vector<uint>> mapModelNodesIds) const;
            void writeTreeWithCorrespondingModels(PhyloTree tree, std::map<uint, vector<uint>> &modelAndNodes) const;
            bool printOutputFileJointLikelihood(const string &fileName, ChromosomeTraitOptimizer* opt) const;
            void createMapOfNodesAndTrait(std::shared_ptr<PhyloTree> stmTree, std::map<uint, string>** mapOfTraitStates, const Alphabet* alpha) const;
            void setIndependentTraitNodeNames(std::map<uint, std::vector<size_t>> &ancestorsTrait, std::map<uint, string>** mapTraitIndependentAncestorNames, const Alphabet* alpha) const;
            std::shared_ptr<LikelihoodCalculationSingleProcess> setTraitLikModel(SingleProcessPhyloLikelihood* lik, ParametrizablePhyloTree* parTree) const;
            void updateMarginalAncestralProbabilitiesAndStates(MarginalAncestralReconstruction *asr, std::map<uint, VVdouble> &posteriorProbs, std::map<uint, vector<size_t>> &mapOfAncestors, const PhyloTree* tree, vector<shared_ptr<PhyloNode>> &nodes);
            void printMarginalAncestralProbabilities(const string &filePath, const PhyloTree* tree, std::map<uint, VVdouble> &posteriorProbs, std::map<uint, vector<size_t>> &mapOfAncestors, const Alphabet* alpha);
            void delete_directory(const std::experimental::filesystem::path& path);
            static bool compare_directory_names(const std::experimental::filesystem::path& a, const std::experimental::filesystem::path& b);
            void rename_simulations_directories(const std::experimental::filesystem::path& dir);

            



        protected:
            void getStatesAndTimes(SiteSimulationResult* simResult, vector<size_t> &states, vector<double> &times, std::shared_ptr<const ParametrizablePhyloTree> parTree, uint nodeId) const;
            void simulateBMTrait(size_t simNum, size_t &count_failed, std::shared_ptr<BrownianMotionSimulator> bmSim);
            void printAncestralStatesInBMSimulation(const string &ancestorsPath, const std::shared_ptr<PhyloTree> tree, const std::unordered_map<uint, double> &statesInNodes);
            void printContinuousTraitFastaFile(const string &countsPath, std::unordered_map<string, double> &statesInLeaves);
            void printOutputFileJointBMLikelihood(const string &fileName, JointChromosomeBMMng* bmMng, ChromosomeNumberOptimizer* independentChromOpt, std::shared_ptr<BrownianMotionLikelihood> bmIndependentLik, bool &nullRejected) const;
            void printJointMLAncestralChromosomeReconstruction(std::shared_ptr<LikelihoodCalculationSingleProcess> likAncestralRec, uint &rootId, int* inferredRootState,  std::map<uint, string>* mapOfTraitStates, ChromosomeTraitOptimizer* traitOpt, ParametrizablePhyloTree* parTree) const;
            void printContinuousTraitBranchStates(const string &fileName, const std::shared_ptr<PhyloTree> tree, const std::unordered_map<uint, double> branchStates);
            
            template<typename T, typename std::enable_if<
                (std::is_same<T, SingleProcessPhyloLikelihood>::value) || 
                (std::is_same<T, JointTraitChromosomeLikelihood>::value)>::type* = nullptr>
            void getTraitThetas(const T* lik, vector<string> &parameterNames, std::unordered_map<string, double> &thetas) const{
                std::regex rgx("theta\\d+");
                auto params = lik->getSubstitutionModelParameters();
                std::smatch match;
                for (size_t i = 0; i < params.size(); i++){
                    if (std::regex_search(parameterNames[i], match, rgx)) {
                        std::string shortParamName = match.str();
                        thetas[shortParamName] = lik->getParameter(parameterNames[i]).getValue();
                    } 
                }
            }
            void getRootFreqsFromTheta(SingleProcessPhyloLikelihood* lik, Vdouble &rootFreqsVec) const;

            void writeTraitMappingPath(std::shared_ptr<PhyloTree> stmTree, std::map<uint, std::vector<size_t>> &ancestors, const string &path) const;
            void writeTraitMappingForNode(std::shared_ptr<PhyloTree> stmTree, std::shared_ptr<PhyloNode> node, std::map<uint, std::vector<size_t>> &ancestors, ofstream &stream) const;
            void removeDummyNodes(std::map<uint, std::vector<size_t>>* ancestors, std::map<uint, std::vector<size_t>>* tempAncestors, const PhyloTree* tree, uint rootId) const;
            VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc, ChromosomeAlphabet* alphaInitial);
            std::map<std::string, std::map<string, double>> extract_alphabet_states(const string &file_path, int &min, int &max, vector<int> &uniqueStates, uint &numberOfComposite);
            void writeRunningParameters(ofstream &outFile) const;
            VectorSiteContainer* createJointTraitChromosomeVscData(VectorSiteContainer* traitVsc) const;
            
            
            std::shared_ptr<LikelihoodCalculationSingleProcess> setHeterogeneousLikInstance(SingleProcessPhyloLikelihood* likProcess, ParametrizablePhyloTree* parTree, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, VectorSiteContainer* vsc, bool forAncestral = false) const;
            std::shared_ptr<NonHomogeneousSubstitutionProcess> setHeterogeneousModel(std::shared_ptr<ParametrizablePhyloTree> tree, SingleProcessPhyloLikelihood* ntl, ValueRef <Eigen::RowVectorXd> rootFreqs,  std::map<int, vector<pair<uint, int>>> sharedParams) const;
            void rescale_tree(PhyloTree* tree, double chrRange);
            void getMaxParsimonyUpperBound(double* parsimonyScore) const;
            // functions to print the tree with ancestral reconstruction
            void printSimulatedDataAndAncestors(SiteSimulationResult* simResult, string &ancestorsPath, string prefix) const;
            void printSimulatedEvoPath(const string outPath, SiteSimulationResult* simResult, bool &success, std::shared_ptr<const ParametrizablePhyloTree> parTree, size_t maxStateIndex = 0) const;
            static string printTree(const PhyloTree& tree);
            static string nodeToParenthesis(const uint nodeId, const PhyloTree& tree);
            std::map<int, vector <double>> getVectorToSetModelParams(SingleProcessPhyloLikelihood* lik, size_t modelIndex = 1) const;
            double getOriginalTreeLength(string &path) const;
            double calculateFreqs(std::unordered_map<string, double> &thetas, size_t &idx) const;
            void setRateChangeTypeToNull() const;
            std::pair<std::shared_ptr<BrownianMotionLikelihood>, ChromosomeNumberOptimizer*> optimizeNullBMModel(JointChromosomeBMMng* bmMng, std::shared_ptr<const std::unordered_map<string, double>> traitData);
            void correctModelParameters(std::map<uint, std::pair<int, std::map<int, vector<double>>>> &complexParamsValues) const;




    };
}
#endif // CHROMEVOL_CHROMOSOMENUMBERMNG_H