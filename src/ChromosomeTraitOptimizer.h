//
// File: ChromosomeNumberOptimizer.h
// Created by: Anat Shafir
// Created on: Wednesday September 2 15:05 2020
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
#ifndef CHROMEVOL_CHROMOSOMETRAITOPTIMIZER_H
#define CHROMEVOL_CHROMOSOMETRAITOPTIMIZER_H

//from bpp-core
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/AbstractNumericalDerivative.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>

//from bpp-seq
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>


//from bpp-phyl

#include <Bpp/Phyl/Tree/PhyloTree.h>
#include "Bpp/Phyl/Likelihood/DataFlow/DataFlowNumeric.h"
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include "UndirectedGraph.h"
#include "ChromosomeSubstitutionModel.h"
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include "ChromEvolOptions.h"
#include "LikelihoodUtils.h"
#include <Bpp/Phyl/OptimizationTools.h>
#include "BaseNumberOptimizer.h"
#include "JointTraitChromosomeLikelihood.h"
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/Character/CharacterSubstitutionModel.h>
// From Seqlib:
#include <vector>
#include <map>
#include <utility>
#include <string>
#include <omp.h>
#include <time.h>
#include <regex>

using namespace std;
namespace bpp
{
    class ChromosomeTraitOptimizer:
    public BaseNumberOptimizer,
    public TraitOptimizer{
    // a class which is used for ChromEvol to run the likelihood optimization procedures
    // with different options available in ChromEvol
    public:
    typedef std::pair<SingleProcessPhyloLikelihood*,SingleProcessPhyloLikelihood*> NullLikelihood;




        protected:
            // for optimizing joint likelihood
            vector <JointTraitChromosomeLikelihood*> vectorOfJointLikelohoods_;
            // for optimizing both likelihood independently
            // both likelihoods will be optimized independently, and finally I will get a product of them
            vector <SingleProcessPhyloLikelihood*> vectorOfLikelihoodsTrait_;
            SingleProcessPhyloLikelihood* optimizedChromosomeLikelihood_;
            const PhyloTree* tree_;
            const ChromosomeAlphabet* alphabetChr_;
            const VectorSiteContainer* vscChr_;
            const VectorSiteContainer* vscTrait_;
            vector<unsigned int> numOfPoints_;
            vector<unsigned int> numOfIterations_;
            uint numOfIterationInBetween_;
            double tolerance_;
            vector<double> *rootFrequenciesChr_;
            vector<double> *rootFrequenciesTrait_;
            std::map<uint, vector<int>> fixedParams_;
            mutable std::map<int, std::vector<std::pair<uint, int>>> sharedParams_;
            string traitModel_;
            size_t numberOfStochasticMappings_;
            bool weightedFreqs_;
            int fixedRootTraitState_;
            std::vector <unsigned int> baseNumberCandidates_;
            double parsimonyBound_;
            vector<string> fixedTraitParams_;
            bool weightedTraitRootFreqs_;
            bool fixedTraitRootFreq_;


            //             vector <SingleProcessPhyloLikelihood*> vectorOfLikelohoods_;
            // //vector <Context> vectorOfContexts_;
            // const PhyloTree* tree_;
            // const ChromosomeAlphabet* alphabet_;
            // const VectorSiteContainer* vsc_;
            // //bool optimizeBaseNumber_;
            // vector<unsigned int> numOfPoints_;
            // vector<unsigned int> numOfIterations_;
            // vector<unsigned int> numOfPointsNextRounds_;
            // vector<unsigned int> numOfIterationsNextRounds_;
            // string typeOfOptimizer_;
            // double tolerance_;
            // bool standardOptimization_;
            // int BrentBracketing_;
            // vector <double> probsForMixedOptimization_;
            // std::map<uint, vector<int>> fixedParams_;
            // mutable std::map<int, std::vector<std::pair<uint, int>>> sharedParams_;
            
            
            

        public:
            ChromosomeTraitOptimizer(
                const PhyloTree* tree,
                const ChromosomeAlphabet* alpha_chr,
                const VectorSiteContainer* vsc_chr,
                const VectorSiteContainer* vsc_trait,
                std::map<uint, uint> baseNumberUpperBound,
                string traitModel,
                size_t numberOfStochasticMappings,
                bool weightedFreqs,
                int fixedRootTraitState,
                double parsimonyBound,
                vector<string> &fixedTraitParams,
                bool weightedTraitRootFreqs,
                bool fixedTraitRootFreq):
                BaseNumberOptimizer(vsc_chr, baseNumberUpperBound),
                    //BaseNumberOptimizer(std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(tempLik_->getLikelihoodCalculation()), optimizeBaseNumber, baseNumOptimizationMethod),
                    vectorOfJointLikelohoods_(),
                    vectorOfLikelihoodsTrait_(),
                    optimizedChromosomeLikelihood_(),
                    //vectorOfContexts_(),
                    tree_(tree),
                    alphabetChr_(alpha_chr),
                    vscChr_(vsc_chr),
                    vscTrait_(vsc_trait),
                    numOfPoints_(),
                    numOfIterations_(),
                    numOfIterationInBetween_(),
                    //baseNumOptimizationMethod_(),
                    //baseNumberUpperBound_(baseNumberUpperBound),
                    tolerance_(),
                    rootFrequenciesChr_(0),
                    rootFrequenciesTrait_(0),
                    fixedParams_(),
                    sharedParams_(ChromEvolOptions::sharedParameters_),
                    traitModel_(traitModel),
                    numberOfStochasticMappings_(numberOfStochasticMappings),
                    weightedFreqs_(weightedFreqs),
                    fixedRootTraitState_(fixedRootTraitState),
                    baseNumberCandidates_(),
                    parsimonyBound_(parsimonyBound),
                    fixedTraitParams_(fixedTraitParams),
                    weightedTraitRootFreqs_(weightedTraitRootFreqs),
                    fixedTraitRootFreq_(fixedTraitRootFreq)
            {
                LikelihoodUtils::updateSharedParameters(sharedParams_, 1, 2);
                for (uint i = 3; i <= static_cast<uint>(ChromEvolOptions::numberOfTraitStates_); i++){
                    LikelihoodUtils::updateSharedParameters(sharedParams_, 1, i);
                }
                fixedParams_[1] = ChromEvolOptions::fixedParams_[1];
                for (uint i = 2; i <= static_cast<uint>(ChromEvolOptions::numberOfTraitStates_); i++){
                    fixedParams_[i] = fixedParams_[1];

                }
            }

            ChromosomeTraitOptimizer(const ChromosomeTraitOptimizer& opt):
                BaseNumberOptimizer(opt),
                vectorOfJointLikelohoods_(opt.vectorOfJointLikelohoods_),
                vectorOfLikelihoodsTrait_(),
                optimizedChromosomeLikelihood_(),
                tree_ (opt.tree_),
                alphabetChr_(opt.alphabetChr_),
                vscChr_(opt.vscChr_),
                vscTrait_(opt.vscTrait_),
                numOfPoints_(opt.numOfPoints_),
                numOfIterations_(opt.numOfIterations_),
                numOfIterationInBetween_(opt.numOfIterationInBetween_),
                tolerance_(opt.tolerance_),
                rootFrequenciesChr_(opt.rootFrequenciesChr_),
                rootFrequenciesTrait_(opt.rootFrequenciesTrait_),
                fixedParams_(opt.fixedParams_),
                sharedParams_(opt.sharedParams_),
                traitModel_(opt.traitModel_),
                numberOfStochasticMappings_(opt.numberOfStochasticMappings_),
                weightedFreqs_(opt.weightedFreqs_),
                fixedRootTraitState_(opt.fixedRootTraitState_),
                baseNumberCandidates_(opt.baseNumberCandidates_),
                parsimonyBound_(opt.parsimonyBound_),
                fixedTraitParams_(opt.fixedTraitParams_),
                weightedTraitRootFreqs_(opt.weightedTraitRootFreqs_),
                fixedTraitRootFreq_(opt.fixedTraitRootFreq_)


            {}
            // ChromosomeTraitOptimizer& operator=(const ChromosomeTraitOptimizer& opt){

            //     vectorOfJointLikelohoods_ = opt.vectorOfJointLikelohoods_;
            //     vectorOfIndependentLikelihoods_ = opt.vectorOfIndependentLikelihoods_;
            //     tree_ = opt.tree_;
            //     alphabetChr_ = opt.alphabetChr_;
            //     vscChr_ = opt.vscChr_;
            //     vscTrait_ = opt.vscTrait_;
            //     numOfPoints_ = opt.numOfPoints_;
            //     numOfIterations_ = opt.numOfIterations_;
            //     numOfIterationInBetween_ = opt.numOfIterationInBetween_;
            //     tolerance_ = opt.tolerance_;
            //     rootFrequenciesChr_ = opt.rootFrequenciesChr_;
            //     rootFrequenciesTrait_ = opt.rootFrequenciesTrait_;
            //     fixedParams_ = opt.fixedParams_;
            //     sharedParams_ = opt.sharedParams_;
            //     traitModel_ = opt.traitModel_;
            //     numberOfStochasticMappings_ = opt.numberOfStochasticMappings_;
            //     weightedFreqs_ = opt.weightedFreqs_;
            //     fixedTraitRootFreq_ = opt.fixedTraitRootFreq_;
            //     baseNumberCandidates_ = opt.baseNumberCandidates_;
            //     parsimonyBound_ = opt.parsimonyBound_;
            //     return *this;
            // }
            ChromosomeTraitOptimizer* clone() const { return new ChromosomeTraitOptimizer(*this); }
            virtual ~ChromosomeTraitOptimizer(){
                //delete vector of joint likelihoods
                clearVectorOfLikelihoods(vectorOfJointLikelohoods_, 0);
                clearVectorOfLikelihoods(vectorOfLikelihoodsTrait_, 0);
                
            }
            //init models
                        // std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelComplexParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds
            //void initModels(std::map<uint, std::pair<int, map<int, std::vector<double>>>> modelComplexParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numberOfModels, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds);
        //     //initialize all the optimization specific members
            void initOptimizer(
                vector<unsigned int> numOfPoints,
                vector<unsigned int> numOfIterations,
                string baseNumOptimizationMethod,
                double tolerance,
                uint numOfIterationInBetween=1)
            {
                numOfPoints_ = numOfPoints;
                numOfIterations_ = numOfIterations;
                baseNumOptimizationMethod_ = baseNumOptimizationMethod;
                tolerance_ = tolerance;
                numOfIterationInBetween_ = numOfIterationInBetween;
                

                getBaseNumCandidates(baseNumberCandidates_, baseNumberUpperBound_);
                setBaseNumberOptimizationMethod(baseNumOptimizationMethod);
                

            }
            const NullLikelihood getNullLikelihood() const{
                return std::make_pair(optimizedChromosomeLikelihood_, vectorOfLikelihoodsTrait_[0]);
            }
            const JointTraitChromosomeLikelihood* getJointLikelihood() const{
                return vectorOfJointLikelohoods_[0];
            }
            void initMultipleLikelihoodPoints(std::map<string, double> &traitModelParams, std::map<uint, std::pair<int, std::map<int, vector<double>>>> &modelParams, const PhyloTree* tree, std::map<uint, uint> baseNumberUpperBound, std::vector<double>* rootFreqsTrait, bool ml);
            vector<double> getFrequenciesFromMapOfParams(std::map<string, double> &traitModelParams, bool random= false);
            void initJointLikelihood(std::map<string, double> traitModelParams, std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, const PhyloTree* tree, std::map<uint, uint> baseNumberUpperBound, std::vector<double>* rootFreqsTrait, bool random, double factor, bool ml);
            static std::map <uint, std::vector<uint>> getNodesForEachModel(std::shared_ptr<PhyloTree> expectedMapping, StochasticMapping* stm);
            void optimizeJointLikelihood();
            void optimizeIndependentLikelihood();
            

            double calculateCriticalLRTClassic(size_t degreesOfFreedom){
                // for some reason here I should write 0.95 instead of 0.05
                return RandomTools::qChisq(0.95, degreesOfFreedom);

            }
            double getLikelihoodNull() const{
                return (optimizedChromosomeLikelihood_->getValue() + vectorOfLikelihoodsTrait_[0]->getValue());
            }

            std::map<string, double> getTraitMLParamsIndependentLik();
            std::map<uint, std::pair<int, std::map<int, vector<double>>>> getChromosomeMLParamsIndependent(uint numOfRequiredModels);
            double getLikelihoodAltrnative() const {
                return vectorOfJointLikelohoods_[0]->getValue();
            }
            size_t getNumberOfParametersJointLikelihood();
            size_t getNumberOfParametersNull();
            double getChiSquare(){
                auto nullLik = getLikelihoodNull();
                auto alternativeLik = getLikelihoodAltrnative();
                return 2*(nullLik-alternativeLik);
            }
            int getThetaIndexIfTheta(const std::string& paramName);
            void setParametersToNewTraitModel(std::map<string, double> &traitModelParams, std::shared_ptr<CharacterSubstitutionModel> traitModel, shared_ptr<IntegerFrequencySet> freqs, bool random);
            void initTraitLikelihoods(std::map<string, double> &traitParams);
            bool testClassicLRT(ofstream& outFile){
                size_t degreesOfFreedom = getNumberOfParametersJointLikelihood()-getNumberOfParametersNull();
                double LRT_c = calculateCriticalLRTClassic(degreesOfFreedom);
                double chiSquare = getChiSquare();
                outFile << "#################################################" << std::endl;
                outFile << "LRT test" << std::endl;
                outFile << "#################################################" << std::endl;
                outFile << "Degrees of freedom: " << degreesOfFreedom <<  std::endl;
                outFile << "chi square statistic is: " << chiSquare << std::endl;
                outFile << "Critical threshold value: " << LRT_c << std::endl;
                if(chiSquare > LRT_c){
                    outFile << "Null model is rejected!" << std::endl;
                }else{
                    outFile << "Null model is not rejected!" << std::endl;
                }
                return (chiSquare > LRT_c);
            }
            void setChromosomeIndependentLikelihood(SingleProcessPhyloLikelihood* optimizedChromosomeLikelihood){
                optimizedChromosomeLikelihood_ = optimizedChromosomeLikelihood;
            }
            const std::map<int, std::vector<pair<uint, int>>> getSharedParams() {return sharedParams_;}
            
            double calculateFreqs(vector<double> &thetas, size_t &idx) const;

            protected:
            void initTraitLikelihood(std::map<string, double> &traitParams, bool random);


            void clearVectorOfLikelihoods(vector<SingleProcessPhyloLikelihood*> &vectorOLikelihoods, size_t new_size){
                while(vectorOLikelihoods.size() > new_size){
                //deleteTreeLikAssociatedAttributes(vectorOfLikelohoods_[vectorOfLikelohoods_.size()-1]);
                    SingleProcessPhyloLikelihood* lik_to_del = vectorOLikelihoods.back(); 
                    vectorOLikelihoods.pop_back();
                    if (lik_to_del){
                        //deleteLikObject(lik_to_del);
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
                }
                vectorOLikelihoods.shrink_to_fit();
            }
            
            void clearVectorOfLikelihoods(vector<JointTraitChromosomeLikelihood*> &vectorOLikelihoods, size_t new_size){
                while(vectorOLikelihoods.size() > new_size){
                //deleteTreeLikAssociatedAttributes(vectorOfLikelohoods_[vectorOfLikelohoods_.size()-1]);
                    JointTraitChromosomeLikelihood* lik_to_del = vectorOLikelihoods.back(); 
                    vectorOLikelihoods.pop_back();
                    if (lik_to_del){
                        //deleteLikObject(lik_to_del);
                        delete lik_to_del;

                    }
                }
                vectorOLikelihoods.shrink_to_fit();
            }


            void printLikelihoodVectorValues(std::vector <JointTraitChromosomeLikelihood*> lik_vec, size_t cycle) const{
                std::cout << "The likelihoods at the end of cycle " + std::to_string(cycle)+" are:" << std::endl;
                for (size_t i = 0; i < lik_vec.size(); i++){
                    std::cout << lik_vec[i]->getValue() << std::endl;
                }

            }
            static bool compareJointLikValues(JointTraitChromosomeLikelihood* lik1, JointTraitChromosomeLikelihood* lik2){
                return (lik1->getValue() < lik2->getValue());
            }

            void getTraitNodes(vector<uint> &modelNodes){
                auto nodes = tree_->getAllNodes();
                for (size_t i = 0; i < nodes.size(); i++){
                    auto nodeId = tree_->getNodeIndex(nodes[i]);
                    if (nodeId == tree_->getRootIndex()){
                        continue;
                    }
                    modelNodes.push_back(nodeId);
                }
            }
            void fillTraitParameters(const string &traitModel, vector<string> &parameterNames, vector<string> &traitParamNames) const;
            


            //std::vector<double> getRootFrequencies(SingleProcessPhyloLikelihood* lik) const;
            //const std::map<int, std::vector<pair<uint, int>>> getSharedParams(){return sharedParams_;}

            //void optimize(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds);
            //vector<SingleProcessPhyloLikelihood*> getVectorOfLikelihoods(){return vectorOfLikelohoods_;}
                    
            
            
            //void writeOutputToFile() const;
            //void printRootFrequencies(SingleProcessPhyloLikelihood* lik, ofstream &outFile) const;
            
            
            


        //protected:
            
            

            
            //void fillVectorOfLikelihoods(SingleProcessPhyloLikelihood* lik, uint numOfIterationsFirstCycle,  size_t currPoint, uint reqNumOfPoints, vector <uint> baseNumCandidates, std::map<int, vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>& fixedParams, vector<SingleProcessPhyloLikelihood*> &vectorOfLiklihoods, string* text, std::map<uint, uint> &baseNumberUpperBounds, omp_lock_t* mutex = 0);
            //void optimizeFirstRound2(SingleProcessPhyloLikelihood* prevLik, uint shiftNode, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, double parsimonyBound, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, vector<uint> numOfPointsNextRounds, vector<uint> numOfIterationsNextRounds, vector<SingleProcessPhyloLikelihood*> &vectorOfLiklihoods, string &text, std::map<uint, uint>* baseNumberBounds, omp_lock_t* mutex= 0);
            

            // add base number map as a parmeter. It might be not safe to use the data member itself (add it to all the optimization functions!!!)


            // // function working on the likelihoods vector object
            //void clearVectorOfLikelihoods(size_t new_size);
            //void clearVectorOfLikelihoods(size_t new_size, std::vector<SingleProcessPhyloLikelihood*> &likelihoodsVec);
            //void deleteLikObject(SingleProcessPhyloLikelihood* lik_to_del);
            


            //print functions
            //void printLikParameters(SingleProcessPhyloLikelihood* lik, unsigned int optimized, string* textToPrint, const string path = "none") const;
           
            //void printLikelihoodVectorValues(vector <SingleProcessPhyloLikelihood*> lik_vec, string* text, size_t cycle) const;

            //void setNewModelAttributes(SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, std::map<uint, pair<int, std::map<int, std::vector<double>>>>* modelParams, std::map<uint, std::vector<uint>>* mapModelNodesIds, std::map<uint, uint>* baseNumberBounds);


            //void updateSharedParameters(std::map<int, vector<std::pair<uint, int>>> &sharedParams, uint prevShift, uint numOfShifts) const;



    };
}
#endif // CHROMEVOL_CHROMOSOMETRAITOPTIMIZER_H

