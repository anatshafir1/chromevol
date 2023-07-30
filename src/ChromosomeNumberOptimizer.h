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
#ifndef CHROMEVOL_CHROMOSOMENUMBEROPTIMIZER_H
#define CHROMEVOL_CHROMOSOMENUMBEROPTIMIZER_H

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
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/JointPhyloLikelihood.h>
#include "BaseNumberOptimizer.h"
#include <Bpp/Phyl/OptimizationTools.h>
// From Seqlib:
#include <vector>
#include <map>
#include <utility>
#include <string>
#include <omp.h>
#include <time.h>
using namespace std;
namespace bpp
{
    class ChromosomeNumberOptimizer:
        public BaseNumberOptimizer{
    // a class which is used for ChromEvol to run the likelihood optimization procedures
    // with different options available in ChromEvol

        protected:
            vector <SingleProcessPhyloLikelihood*> vectorOfLikelohoods_;
            JointPhyloLikelihood* singleLikelihood_;
            //vector <Context> vectorOfContexts_;
            const PhyloTree* tree_;
            const ChromosomeAlphabet* alphabet_;
            const VectorSiteContainer* vsc_;
            //bool optimizeBaseNumber_;
            vector<unsigned int> numOfPoints_;
            vector<unsigned int> numOfIterations_;
            vector<unsigned int> numOfPointsNextRounds_;
            vector<unsigned int> numOfIterationsNextRounds_;
            string typeOfOptimizer_;
            //string baseNumOptimizationMethod_;
            //mutable std::map<uint, uint> baseNumberUpperBound_;
            double tolerance_;
            bool standardOptimization_;
            int BrentBracketing_;
            vector <double> probsForMixedOptimization_;
            std::map<uint, vector<int>> fixedParams_;
            mutable std::map<int, std::vector<std::pair<uint, int>>> sharedParams_;
            uint numOfShiftsForward_;
            bool backwardPhaseStarted_;
            // number of shifts. For each shift, a new added node
            std::map<uint, vector<uint>> prevModelsPartitions_;
            std::map<uint, std::vector<std::pair<string, double>>> prevModelParams_;
            std::map<uint, std::pair<double, double>> prevModelsAICcLikValues_;
            std::map<uint, std::vector<double>> prevModelsRootFrequencies_;
            
            
            

        public:
            ChromosomeNumberOptimizer(
                const PhyloTree* tree,
                const ChromosomeAlphabet* alpha,
                const VectorSiteContainer* vsc,
                std::map<uint, uint> baseNumberUpperBound):
                    BaseNumberOptimizer(vsc, baseNumberUpperBound),
                    vectorOfLikelohoods_(),
                    singleLikelihood_(0),
                    //vectorOfContexts_(),
                    tree_(tree),
                    alphabet_(alpha),
                    vsc_(vsc),
                    numOfPoints_(),
                    numOfIterations_(),
                    numOfPointsNextRounds_(),
                    numOfIterationsNextRounds_(),
                    typeOfOptimizer_(),  
                    tolerance_(),
                    standardOptimization_(),
                    BrentBracketing_(),
                    probsForMixedOptimization_(),
                    fixedParams_(),
                    sharedParams_(),
                    numOfShiftsForward_(),
                    backwardPhaseStarted_(),
                    prevModelsPartitions_(),
                    prevModelParams_(),
                    prevModelsAICcLikValues_(),
                    prevModelsRootFrequencies_()
            {}
            // when using an optimizer for joint likelihood
            ChromosomeNumberOptimizer(
                const JointPhyloLikelihood* jointLik,
                const ChromosomeAlphabet* alpha,
                const VectorSiteContainer* vsc,
                std::map<uint, uint> &baseNumberUpperBound,
                std::map<int, std::vector<std::pair<uint, int>>> &sharedParams,
                std::map<uint, vector<int>> &fixedParams):
                    BaseNumberOptimizer(vsc, baseNumberUpperBound),
                    vectorOfLikelohoods_(),
                    singleLikelihood_(0),
                    //vectorOfContexts_(),
                    tree_(0), // do we need a tree here
                    alphabet_(alpha),
                    vsc_(0), // is it needed?
                    numOfPoints_(),
                    numOfIterations_(),
                    numOfPointsNextRounds_(),
                    numOfIterationsNextRounds_(),
                    typeOfOptimizer_(),  
                    tolerance_(),
                    standardOptimization_(),
                    BrentBracketing_(),
                    probsForMixedOptimization_(),
                    fixedParams_(fixedParams),
                    sharedParams_(sharedParams),
                    numOfShiftsForward_(),
                    backwardPhaseStarted_(),
                    prevModelsPartitions_(),
                    prevModelParams_(),
                    prevModelsAICcLikValues_(),
                    prevModelsRootFrequencies_()
            {}

            ChromosomeNumberOptimizer(const ChromosomeNumberOptimizer& opt):
                BaseNumberOptimizer(opt),
                vectorOfLikelohoods_(opt.vectorOfLikelohoods_),
                singleLikelihood_(opt.singleLikelihood_),
                //vectorOfContexts_(opt.vectorOfContexts_),
                tree_ (opt.tree_),
                alphabet_(opt.alphabet_),
                vsc_(opt.vsc_),
                numOfPoints_(opt.numOfPoints_),
                numOfIterations_(opt.numOfIterations_),
                numOfPointsNextRounds_(opt.numOfPointsNextRounds_),
                numOfIterationsNextRounds_(opt.numOfIterationsNextRounds_),
                typeOfOptimizer_(opt.typeOfOptimizer_),
                tolerance_(opt.tolerance_),
                standardOptimization_(opt.standardOptimization_),
                BrentBracketing_(opt.BrentBracketing_),
                probsForMixedOptimization_(opt.probsForMixedOptimization_),
                fixedParams_(opt.fixedParams_),
                sharedParams_(opt.sharedParams_),
                numOfShiftsForward_(opt.numOfShiftsForward_),
                backwardPhaseStarted_(opt.backwardPhaseStarted_),
                prevModelsPartitions_(opt.prevModelsPartitions_),
                prevModelParams_(opt.prevModelParams_),
                prevModelsAICcLikValues_(opt.prevModelsAICcLikValues_),
                prevModelsRootFrequencies_(opt.prevModelsRootFrequencies_)

            {}
            ChromosomeNumberOptimizer& operator=(const ChromosomeNumberOptimizer& opt){
                vectorOfLikelohoods_ = opt.vectorOfLikelohoods_;
                singleLikelihood_ = opt.singleLikelihood_;
                //vectorOfContexts_ = opt.vectorOfContexts_;
                tree_ = opt.tree_;
                alphabet_ = opt.alphabet_;
                vsc_ = opt.vsc_;
                numOfPoints_ = opt.numOfPoints_;
                numOfIterations_ = opt.numOfIterations_;
                numOfPointsNextRounds_ = opt.numOfPointsNextRounds_;
                numOfIterationsNextRounds_ = opt.numOfIterationsNextRounds_;
                typeOfOptimizer_ = opt.typeOfOptimizer_;
                tolerance_ = opt.tolerance_;
                standardOptimization_ = opt.standardOptimization_;
                BrentBracketing_ = opt.BrentBracketing_;
                probsForMixedOptimization_ = opt.probsForMixedOptimization_;
                fixedParams_ = opt.fixedParams_;
                sharedParams_ = opt.sharedParams_;
                numOfShiftsForward_ = opt.numOfShiftsForward_;
                backwardPhaseStarted_ = opt.backwardPhaseStarted_;
                prevModelsPartitions_ = opt.prevModelsPartitions_;
                prevModelParams_ = opt.prevModelParams_;
                prevModelsAICcLikValues_ = opt.prevModelsAICcLikValues_;
                prevModelsRootFrequencies_ = opt.prevModelsRootFrequencies_;
                return *this;
            }
            ChromosomeNumberOptimizer* clone() const { return new ChromosomeNumberOptimizer(*this); }
            virtual ~ChromosomeNumberOptimizer(){clearVectorOfLikelihoods(0);};
            //init models
                        // std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelComplexParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds
            //void initModels(std::map<uint, std::pair<int, map<int, std::vector<double>>>> modelComplexParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numberOfModels, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds);
        //     //initialize all the optimization specific members
            void initOptimizer(
                vector<unsigned int> numOfPoints,
                vector<unsigned int> numOfIterations,
                vector<unsigned int> numOfPointsNextRounds,
                vector<unsigned int> numOfIterationsNextRounds,
                string typeOfOptimizer,
                string baseNumOptimizationMethod,
                double tolerance,
                bool standardOptimization,
                int BrentBracketing,
                vector <double>& probsForMixedOptimization)
            {
                numOfPoints_ = numOfPoints;
                numOfIterations_ = numOfIterations;
                numOfPointsNextRounds_ = numOfPointsNextRounds;
                numOfIterationsNextRounds_ = numOfIterationsNextRounds;
                typeOfOptimizer_ = typeOfOptimizer;
                tolerance_ = tolerance;
                standardOptimization_ = standardOptimization;
                BrentBracketing_ =BrentBracketing;
                probsForMixedOptimization_ = probsForMixedOptimization;
                backwardPhaseStarted_ = false;
                setBaseNumberOptimizationMethod(baseNumOptimizationMethod);
                

            }
            // typename std::enable_if<(std::is_same<U, Eigen::MatrixXd>::value) && (std::is_same<V, Eigen::MatrixXd>::value), void>::type
            //template<typename T>
            //typename std::enable_if<(std::is_same<T, SingleProcessPhyloLikelihood>::value) || (std::is_same<T, JointPhyloLikelihood>::value),unsigned int>::type
            //optimizeModelParameters(T* tl, double tol, unsigned int maxNumOfIterations, vector<unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, std::vector<string>* chromosomeParamNames = 0);
            const std::map<uint, vector<uint>> getPreviousModelsPartitions() const{return prevModelsPartitions_;}
            const std::map<uint, std::vector<std::pair<string, double>>> getPreviousModelsParameters() const{return prevModelParams_;}
            const std::map<uint, std::pair<double, double>> getPreviousModelsAICcValues() const{return prevModelsAICcLikValues_;}
            const std::map <uint,std::vector<double>> getPrevModelsRootFreqs() const{return prevModelsRootFrequencies_;}
            std::vector<double> getRootFrequencies(SingleProcessPhyloLikelihood* lik) const;
            void getParameterNamesAndValues(SingleProcessPhyloLikelihood* lik, uint numOfModels);
            void setInitialModelRepresentitives(std::map<uint, vector<uint>> &initialPartition);
            const std::map<int, std::vector<pair<uint, int>>> getSharedParams(){return sharedParams_;}
            const double getAICOfBestModel() const {
                std::map<uint, vector<int>> fixedParams = fixedParams_;
                size_t numOfFixedParams = LikelihoodUtils::getNumberOfFixedParams(vectorOfLikelohoods_[0], fixedParams);
                return calculateModelSelectionCriterion(vectorOfLikelohoods_[0], numOfFixedParams);
            }
            void runNewBranchModel(omp_lock_t &mutex, SingleProcessPhyloLikelihood* lik, vector<uint> &candidateShiftNodesIds, size_t i, uint numOfShifts, double parsimonyBound, uint numOfPoints, SingleProcessPhyloLikelihood** bestCandidateLik, double* bestAICc, uint* minAICcNode);
            void optimizeMultiProcessModel(std::map<int, std::vector<pair<uint, int>>>* sharedParams ,std::map<uint, vector<int>>* fixedParams, vector<unsigned int> &numOfPoints, vector<unsigned int> &numOfIterations, std::map<uint, uint> &baseNumberUpperBounds, vector<SingleProcessPhyloLikelihood*>* perCandidateLikVec, string* textToPrint, omp_lock_t* mutex = 0);
            void optimize(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds);
            void optimizeInParallel(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds);
            void optimizeBackwards(double maxParsimony, bool parallel);
            vector<SingleProcessPhyloLikelihood*> getVectorOfLikelihoods(){return vectorOfLikelohoods_;}
                    
            static void getMapOfModelAndNodeIdsBackward(std::map<uint, vector<uint>> &mapModelNodesIds, SingleProcessPhyloLikelihood* lik, std::map<uint, uint> &modelsMap, std::map<uint, vector<uint>> &modelForMerge, uint rootId = 0);
            
            
            //void writeOutputToFile() const;
            void printRootFrequencies(SingleProcessPhyloLikelihood* lik, ofstream &outFile) const;
            
            
            
            
            void setIterNumForNextRound(std::vector<uint> iterNum){
                numOfIterationsNextRounds_ = iterNum;
            }
            void setPointsNumForNextRound(std::vector<uint> pointsNum){
                numOfPointsNextRounds_ = pointsNum;
            }




        protected:
            
            
        //     // //functions of optimization
            void printLog(string* fileOut, string text){
                if (fileOut){
                    *fileOut += (text);
                }else{
                    std::cout << text;
                }
            }
            
            void fillVectorOfLikelihoods(SingleProcessPhyloLikelihood* lik, uint numOfIterationsFirstCycle,  size_t currPoint, uint reqNumOfPoints, vector <uint> baseNumCandidates, std::map<int, vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>& fixedParams, vector<SingleProcessPhyloLikelihood*> &vectorOfLiklihoods, string* text, std::map<uint, uint> &baseNumberUpperBounds, omp_lock_t* mutex = 0);
            void initLikelihoods(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds, uint numOfModels, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams);
            void optimizeFirstRound(std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, double parsimonyBound, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, vector<uint> numOfPointsNextRounds, vector<uint> numOfIterationsNextRounds, vector<SingleProcessPhyloLikelihood*> &vectorOfLiklihoods, string* text, std::map<uint, uint>* baseNumberBounds, std::map<uint, uint>* mapOfModelsBackward, std::map<uint, pair<int, std::map<int, std::vector<double>>>>* prevModelParamsBackward, std::map<uint, vector<uint>>* modelsBackwards, omp_lock_t* mutex = 0);
            //void optimizeFirstRound2(SingleProcessPhyloLikelihood* prevLik, uint shiftNode, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, double parsimonyBound, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, vector<uint> numOfPointsNextRounds, vector<uint> numOfIterationsNextRounds, vector<SingleProcessPhyloLikelihood*> &vectorOfLiklihoods, string &text, std::map<uint, uint>* baseNumberBounds, omp_lock_t* mutex= 0);
            

            // add base number map as a parmeter. It might be not safe to use the data member itself (add it to all the optimization functions!!!)
            
            //template<typename T> unsigned int optimizeModelParametersOneDimension(T* tl, double tol, unsigned int maxNumOfIterations, std::vector<unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, std::vector<string>* chromosomeParamNames, bool mixed = false, unsigned int currentIterNum = 0);
            //template<typename T> unsigned int optimizeMultiDimensions(T* tl, double tol, unsigned int maxNumOfIterations, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, std::vector<string>* chromosomeParamNames, bool mixed = false, unsigned int currentIterNum = 0);

            //template<typename T> unsigned int useMixedOptimizers(T* tl, double tol, unsigned int maxNumOfIterations, vector <unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, std::vector<string>* chromosomeParamNames);

            // // function working on the likelihoods vector object
            void clearVectorOfLikelihoods(size_t new_size);
            void clearVectorOfLikelihoods(size_t new_size, std::vector<SingleProcessPhyloLikelihood*> &likelihoodsVec);
            void deleteLikObject(SingleProcessPhyloLikelihood* lik_to_del);
            

            // // helper functions for optimization
            void checkLegalUseOfGradientOptimization();

            
            //vector <string> getNonFixedParams(SingleProcessPhyloLikelihood* tl, ParameterList &allParams, map<uint, vector<int>>* fixedParams, vector<string>* chromosomeParamNames) const;
            //string findParameterNameInModel(string fullParameterName) const;
            //void setNewBounds(const ParameterList params, Parameter &param, map<string, pair<string, bool>> &paramPairsMap, double* lowerBound, const ChromosomeSubstitutionModel* model);

            //print functions
            //template<typename T> void printLikParameters(T* lik, unsigned int optimized, string* textToPrint, const string path = "none") const;
           
            void printLikelihoodVectorValues(vector <SingleProcessPhyloLikelihood*> lik_vec, string* text, size_t cycle) const;

            /*********************************************************
             * Functions associated with heterogenous ChromEvol models
            **********************************************************/
            double calculateModelSelectionCriterion(SingleProcessPhyloLikelihood* lik, size_t numOfFixedParams) const;
            double calculateAICc(SingleProcessPhyloLikelihood* lik, size_t numOfFixedParams) const;
            double calculateAIC(SingleProcessPhyloLikelihood* lik, size_t numOfFixedParams) const;
            //getNewLikObject(SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<int, uint>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, uint numOfPoints, std::map<uint, vector<int>> &fixedParams, double parsimonyBound)
            void setNewModelAttributes(SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, std::map<uint, pair<int, std::map<int, std::vector<double>>>>* modelParams, std::map<uint, std::vector<uint>>* mapModelNodesIds, std::map<uint, uint>* baseNumberBounds);

            SingleProcessPhyloLikelihood* getSingleNewLikObject(std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, double parsimonyBound, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, uint iteration, std::map<uint, uint>* baseNumberBounds);
            void getNewLikObjectForParallelRuns(std::vector<SingleProcessPhyloLikelihood*> &perCandidateLikVec, SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, uint numOfPoints, std::map<uint, vector<int>> &fixedParams, double parsimonyBound);
            void getNewLikObject(SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, uint numOfPoints, std::map<uint, vector<int>> &fixedParams, double parsimonyBound);

            void ifNanTryToResampleLikObject(SingleProcessPhyloLikelihood** lik, const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, double parsimonyBound, int numOfPoints, std::map<uint, vector<int>> &fixedParams, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams);
            //void optimizeSingleHeterogeneousModel(size_t index, int maxNumOfModels, std::vector<uint> &candidateShiftNodesIds, vector<uint> &baseNumCandidates);
            void getValidCandidatesForShift(std::vector<uint> &candidateShiftNodesIds, int minCladeSize, uint numOfShifts);
            size_t getValidCandidatesForShiftRec(uint nodeId, std::vector<uint> &candidateShiftNodesIds, int minCladeSize, vector<uint> &shifting_nodes, bool shifting_node);
            

            
            /*********************************************************
             * Functions for the Backward phase
            **********************************************************/
          void fillWithFathers(vector<uint> &fathers, vector<uint> &nodes);
          void optimizeMergedModels(SingleProcessPhyloLikelihood* finalLikBackward, std::map<uint, vector<uint>> &modelsToBeMerged, std::map<std::pair<uint, uint>, double>* pairsOfLikelihoods, double maxParsimony, omp_lock_t* mutex);
          size_t getMaxNumOfMergingModels(std::map<uint, vector<uint>> &modelsToMerge);
          void mergeModels(std::map<uint, vector<uint>> modelsToMerge, SingleProcessPhyloLikelihood* lik, std::map<uint, vector<int>> &fixedParams, std::map<int, std::vector<std::pair<uint, int>>> &updatedSharedParams, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParamsPrevModel, std::map<uint, uint> &modelNums, std::map<uint, uint> &baseNumberBounds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams);
          void updateSharedParametersBackwards(std::map<int, std::vector<std::pair<uint, int>>> &sharedParams, std::map<int, std::vector<std::pair<uint, int>>> &updatedSharedParams, std::map<uint,uint> &mapOfModels);
          SingleProcessPhyloLikelihood* getBackwardLikObject(std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, std::map<uint, vector<int>> &fixedParams, double parsimonyBound, std::map<uint, std::vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &prevModelParamsBackward, std::map<uint, uint> &mapOfModelsBackward, std::map<uint, vector<uint>> &models, uint iteration, std::map<uint, uint>* baseNumberBounds);
          pair<int, std::map<int, std::vector<double>>> getMeanParameters(std::map<uint, pair<int, std::map<int, std::vector<double>>>> &prevModelParamsBackward, vector<uint> modelsBackwards, uint modelForMerge);
          
          void getMapOfMergedModels(std::map<uint, uint> &mapOfModels, std::map<uint, vector<uint>> &modelsToMerge, SingleProcessPhyloLikelihood* prevLik);
          void mergeMultipleModelClusters(SingleProcessPhyloLikelihood* finalLikBackward, std::map<uint, vector<uint>> &rootAndVerticesToMerge, double maxParsimony);
          //mergeMultipleModelClusters(SingleProcessPhyloLikelihood* prevLik, std::map<uint, vector<uint>> &rootAndVerticesToMerge, double maxParsimony)
          bool isModelADirectSubtreeOfAnother(SingleProcessPhyloLikelihood* lik, size_t indexModel1, size_t indexModel2);
          /**************************************************************************************************************************************/

          // functions related to both joint model and single likelihood objects

            template <typename T>
            typename std::enable_if<(std::is_same<T, SingleProcessPhyloLikelihood>::value) || (std::is_same<T, JointPhyloLikelihood>::value),void>::type
            printLikParameters(T* lik, unsigned int optimized, string* textToPrint, const string filePath = "none") const{
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
                std::vector<std::string> paramsNames;
                auto substitutionModelParams = lik->getSubstitutionModelParameters();
                if (tree_){
                    paramsNames = substitutionModelParams.getParameterNames();

                }else{
                    LikelihoodUtils::fixSuffixForJointLikParamNames(substitutionModelParams, paramsNames);

                }
                
                for (int i = 0; i < (int)(paramsNames.size()); i++){
                    if (paramsNames[i].find("Chromosome.baseNum_") != std::string::npos){
                        //text += paramsNames[i]+ " = "+ std::to_string((int)(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()))+"\n";
                        text += paramsNames[i]+ " = "+ std::to_string((int)(lik->getParameters().getParameter(paramsNames[i]).getValue()))+"\n";
                        if (filePath != "none"){
                            //outFile <<  paramsNames[i] << " = "<< (int)(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()) <<endl;
                            outFile <<  paramsNames[i] << " = "<< (int)(lik->getParameters().getParameter(paramsNames[i]).getValue()) <<endl;
                        }
                    }else{
                        //text += (paramsNames[i] + " = "+ std::to_string(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()) + "\n");
                        text += (paramsNames[i] + " = "+ std::to_string(lik->getParameters().getParameter(paramsNames[i]).getValue()) + "\n");
                        if (filePath != "none"){
                            //outFile << paramsNames[i] << " = "<< lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue() <<endl;
                            outFile << paramsNames[i] << " = "<< lik->getParameters().getParameter(paramsNames[i]).getValue() <<endl;
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
            template <typename T>
            typename std::enable_if<(std::is_same<T, SingleProcessPhyloLikelihood>::value), uint>::type
            getNumberOfModelsInLikObject(T* likObject) const{
                return static_cast<uint>(likObject->getSubstitutionProcess().getNumberOfModels());
            }
            template <typename T>
            typename std::enable_if<(std::is_same<T, JointPhyloLikelihood>::value), uint>::type
            getNumberOfModelsInLikObject(T* likObject) const{
                return static_cast<uint>(likObject->getPhylo2()->getSubstitutionProcess().getNumberOfModels());
            }

            template <typename T>
            typename std::enable_if<(std::is_same<T, SingleProcessPhyloLikelihood>::value) || (std::is_same<T, JointPhyloLikelihood>::value),vector<string>>::type
            getNonFixedParams(T* tl, ParameterList &allParams, std::map<uint, vector<int>>* fixedParams, vector<string>* chromosomeParamNames) const{
                std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
                std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model
                if (chromosomeParamNames){
                    LikelihoodUtils::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, *chromosomeParamNames, &sharedParams_);
                }else{
                    
                    LikelihoodUtils::updateMapsOfParamTypesAndNames<T>(typeWithParamNames, &paramNameAndType, tl);
                }
                uint numOfModels = getNumberOfModelsInLikObject<T>(tl);
                // if constexpr (std::is_same_v<T, SingleProcessPhyloLikelihood>){
                //     numOfModels = static_cast<uint>(tl->getSubstitutionProcess().getNumberOfModels());

                // }else if constexpr (std::is_same_v<T, JointPhyloLikelihood>){
                //     numOfModels = static_cast<uint>(tl->getPhylo2()->getSubstitutionProcess().getNumberOfModels());

                // }else{
                //     throw Exception("ERROR!!!! ChromosomeNumberOptimizer::getNonFixedParams(): Type that substitute T is not correct!!!");
                // }
                
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

        public:
            template <typename T>
            typename std::enable_if<(std::is_same<T, SingleProcessPhyloLikelihood>::value) || (std::is_same<T, JointPhyloLikelihood>::value),unsigned int>::type
            optimizeMultiDimensions(T* tl, double tol, unsigned int maxNumOfIterations, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, std::vector<string>* chromosomeParamNames, bool mixed=false, unsigned int currentIterNum=0){
                DerivableSecondOrder* f = tl;
                unique_ptr<AbstractNumericalDerivative> fnum;
                fnum.reset(new TwoPointsNumericalDerivative(f));
                fnum->setInterval(0.0000001);
                ConjugateGradientMultiDimensions* optimizer = new ConjugateGradientMultiDimensions(fnum.get());
                if (chromosomeParamNames){
                    fnum->setParametersToDerivate(*chromosomeParamNames);
                }else{
                    fnum->setParametersToDerivate(tl->getSubstitutionModelParameters().getParameterNames());
                }
                optimizer->setVerbose(1);
                optimizer->setProfiler(0);
                optimizer->setMessageHandler(0);
                optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
                optimizer->getStopCondition()->setTolerance(tol* 0.1);
                optimizer->setMaximumNumberOfEvaluations(1000);
                std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
                std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model
                if (chromosomeParamNames){
                    LikelihoodUtils::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, *chromosomeParamNames, &sharedParams_);
                }else{
                    LikelihoodUtils::updateMapsOfParamTypesAndNames<T>(typeWithParamNames, &paramNameAndType, tl);
                }
                size_t startCompositeParams = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();

                string text = "";
                unsigned int numOfEvaluations = 0;
                double currentLikelihood = tl->getValue();
                double prevLikelihood;
                for (size_t i = 0; i < maxNumOfIterations; i++){
                    if (!chromosomeParamNames){
                        if(mixed){
                            text += ("Iteration #"+ std::to_string(currentIterNum)+"\n");

                        }else{
                            text += ("Iteration #"+ std::to_string(i)+"\n");
                        }

                    }

        
                    ParameterList paramsFull = tl->getSubstitutionModelParameters();
                    std::vector <string> nonFixedparamsNames = getNonFixedParams(tl, paramsFull, fixedParams, chromosomeParamNames);
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
                    printLikParameters<T>(tl, 1, textToPrint);
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



            /********************************************************************************************************************************************************/
            template <typename T>
            typename std::enable_if<(std::is_same<T, SingleProcessPhyloLikelihood>::value) || (std::is_same<T, JointPhyloLikelihood>::value),unsigned int>::type
            optimizeModelParametersOneDimension(T* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates, std::map<int, vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberBounds, std::vector<string>* chromosomeParamNames, bool mixed=false, unsigned curentIterNum=0){

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
            // updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>>* paramNameAndType, std::vector<std::string> namesAllParams, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::string suffix)
            if (!(chromosomeParamNames)){
                LikelihoodUtils::updateMapsOfParamTypesAndNames<T>(typeWithParamNames, &paramNameAndType, tl);

            }else{
                LikelihoodUtils::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, *chromosomeParamNames, &sharedParams_);

            }

            ParameterList params;
            vector<string> parametersNames;
            if (chromosomeParamNames){
                parametersNames = *chromosomeParamNames;
            }else{
                parametersNames = tl->getSubstitutionModelParameters().getParameterNames();
            }
            size_t nbParams = parametersNames.size();
            // starting iterations of optimization
            for (size_t i = 0; i < maxNumOfIterations; i++){
                if (!chromosomeParamNames){
                    if (mixed){
                        text = "Iteration #" + std::to_string(curentIterNum)+ "\n";
        

                    }else{
                        text = "Iteration #" +std::to_string(i) + "\n";
                    }

                }

                printLog(textToPrint, text);
                //ParameterList params = tl->getParameters();// = tl->getParameters();
                //ParameterList substitutionModelParams = tl->getSubstitutionModelParameters();
    
                prevLikelihood = currentLikelihood;
    
                for (size_t j = 0; j < nbParams; j ++){
                    params = tl->getParameters();
                    const string nameOfParam = parametersNames[j];
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
                        std::shared_ptr<IntervalConstraint> intervalFuncUpdatedTL = dynamic_pointer_cast<IntervalConstraint>(tl->getParameter(nameOfParam).getConstraint());

                    }else{
                        // baseNumber parameter
                        if (baseNumOptimizationMethod_ != "Brent"){
                            if (!std::count((*fixedParams)[paramNameAndType[nameOfParam].second].begin(), (*fixedParams)[paramNameAndType[nameOfParam].second].end(), ChromosomeSubstitutionModel::BASENUM)){
                                optimizeBaseNum(tl, parametersNames, j, baseNumCandidates, &currentLikelihood, lowerBound, upperBound, nameOfParam, params, paramNameAndType[nameOfParam].second, baseNumberBounds);
                                //text = "parameter value after optimization "+ std::to_string(tl->getLikelihoodCalculation()->getParameter(param.getName()).getValue())+ "\n";
                                text = "parameter value after optimization "+ std::to_string(tl->getParameters().getParameter(param.getName()).getValue())+ "\n";
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
                    //text = "parameter value after optimization "+ std::to_string(tl->getLikelihoodCalculation()->getParameter(param.getName()).getValue())+ "\n";
                    text = "parameter value after optimization "+ std::to_string(tl->getParameters().getParameter(param.getName()).getValue())+ "\n";
                    printLog(textToPrint, text);
                    text = "***\n";
                    printLog(textToPrint, text);
                }
                printLikParameters<T>(tl, 1, textToPrint);
    
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
        /**********************************************************************************************************************************************************/
        template <typename T>
        typename std::enable_if<(std::is_same<T, SingleProcessPhyloLikelihood>::value) || (std::is_same<T, JointPhyloLikelihood>::value),unsigned int>::type
        useMixedOptimizers(T* tl, double tol, unsigned int maxNumOfIterations, vector <unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, std::vector<string>* chromosomeParamNames){
            std::vector<size_t> optimization = RandomTools::randMultinomial(maxNumOfIterations, probsForMixedOptimization_);
            unsigned int numOfEvaluations = 0;
            for (size_t i = 0; i < maxNumOfIterations; i++){
                double prevLikelihood = tl->getValue();
                if (optimization[i] == 0){
                    string text = "Optimizing with Brent\n";
                    printLog(textToPrint, text);
                    numOfEvaluations += optimizeModelParametersOneDimension<T>(tl, tol, 1, baseNumCandidates, sharedParams, fixedParams, textToPrint, baseNumberUpperBounds, chromosomeParamNames, true, (unsigned int)i);
                }else{
                    string text = "Optimizing with Gradient Descent\n";
                    printLog(textToPrint, text);
                    numOfEvaluations += optimizeMultiDimensions<T>(tl, tol, 1, sharedParams, fixedParams, textToPrint, baseNumberUpperBounds, chromosomeParamNames, true, (unsigned int)i);
                }
                double currentLikValue = tl->getValue();
                if (std::abs(prevLikelihood-currentLikValue) < tol){
                    break;
                }


            }
    
            return numOfEvaluations;

        }

        /**********************************************************************************************************************************************************/


        template <typename T>
        typename std::enable_if<(std::is_same<T, SingleProcessPhyloLikelihood>::value) || (std::is_same<T, JointPhyloLikelihood>::value),unsigned int>::type
        optimizeModelParameters(T* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, std::vector<string>* chromosomeParamNames=0){
            unsigned int numOfEvaluations = 0;
 
            if (typeOfOptimizer_ == "Brent"){
                numOfEvaluations += optimizeModelParametersOneDimension<T>(tl, tol, maxNumOfIterations, baseNumCandidates, sharedParams, fixedParams, textToPrint, baseNumberUpperBounds, chromosomeParamNames);
            }else if (typeOfOptimizer_ == "gradient"){
                checkLegalUseOfGradientOptimization();
                numOfEvaluations += optimizeMultiDimensions<T>(tl, tol, maxNumOfIterations, sharedParams, fixedParams, textToPrint, baseNumberUpperBounds, chromosomeParamNames);

            }else{
                checkLegalUseOfGradientOptimization();
                numOfEvaluations += useMixedOptimizers<T>(tl, tol, maxNumOfIterations, baseNumCandidates, sharedParams, fixedParams, textToPrint, baseNumberUpperBounds, chromosomeParamNames);
            }
        
            return numOfEvaluations;
    
        }

    };
}
#endif // CHROMEVOL_CHROMOSOMENUMBEROPTIMIZER_H