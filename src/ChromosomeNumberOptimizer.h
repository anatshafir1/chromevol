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
    class ChromosomeNumberOptimizer{
    // a class which is used for ChromEvol to run the likelihood optimization procedures
    // with different options available in ChromEvol

        private:
            vector <SingleProcessPhyloLikelihood*> vectorOfLikelohoods_;
            //vector <Context> vectorOfContexts_;
            const PhyloTree* tree_;
            const ChromosomeAlphabet* alphabet_;
            const VectorSiteContainer* vsc_;
            bool optimizeBaseNumber_;
            vector<unsigned int> numOfPoints_;
            vector<unsigned int> numOfIterations_;
            vector<unsigned int> numOfPointsNextRounds_;
            vector<unsigned int> numOfIterationsNextRounds_;
            string typeOfOptimizer_;
            string baseNumOptimizationMethod_;
            mutable std::map<uint, uint> baseNumberUpperBound_;
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
                    vectorOfLikelohoods_(),
                    //vectorOfContexts_(),
                    tree_(tree),
                    alphabet_(alpha),
                    vsc_(vsc),
                    optimizeBaseNumber_(),
                    numOfPoints_(),
                    numOfIterations_(),
                    numOfPointsNextRounds_(),
                    numOfIterationsNextRounds_(),
                    typeOfOptimizer_(),
                    baseNumOptimizationMethod_(),
                    baseNumberUpperBound_(baseNumberUpperBound),
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

            ChromosomeNumberOptimizer(const ChromosomeNumberOptimizer& opt):
                vectorOfLikelohoods_(opt.vectorOfLikelohoods_),
                //vectorOfContexts_(opt.vectorOfContexts_),
                tree_ (opt.tree_),
                alphabet_(opt.alphabet_),
                vsc_(opt.vsc_),
                optimizeBaseNumber_(opt.optimizeBaseNumber_),
                numOfPoints_(opt.numOfPoints_),
                numOfIterations_(opt.numOfIterations_),
                numOfPointsNextRounds_(opt.numOfPointsNextRounds_),
                numOfIterationsNextRounds_(opt.numOfIterationsNextRounds_),
                typeOfOptimizer_(opt.typeOfOptimizer_),
                baseNumOptimizationMethod_(opt.baseNumOptimizationMethod_),
                baseNumberUpperBound_(opt.baseNumberUpperBound_),
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
                //vectorOfContexts_ = opt.vectorOfContexts_;
                tree_ = opt.tree_;
                alphabet_ = opt.alphabet_;
                vsc_ = opt.vsc_;
                optimizeBaseNumber_ = opt.optimizeBaseNumber_;
                numOfPoints_ = opt.numOfPoints_;
                numOfIterations_ = opt.numOfIterations_;
                numOfPointsNextRounds_ = opt.numOfPointsNextRounds_;
                numOfIterationsNextRounds_ = opt.numOfIterationsNextRounds_;
                typeOfOptimizer_ = opt.typeOfOptimizer_;
                baseNumOptimizationMethod_ = opt.baseNumOptimizationMethod_;
                baseNumberUpperBound_ = opt.baseNumberUpperBound_;
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
                baseNumOptimizationMethod_ = baseNumOptimizationMethod;
                tolerance_ = tolerance;
                standardOptimization_ = standardOptimization;
                BrentBracketing_ =BrentBracketing;
                probsForMixedOptimization_ = probsForMixedOptimization;
                backwardPhaseStarted_ = false;
                

            }
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
            unsigned int optimizeModelParameters(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, vector<unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds);//, unsigned int inwardBracketing, bool standardOptimization);
            unsigned int optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector<unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, bool mixed = false, unsigned int currentIterNum = 0);
            unsigned int optimizeMultiDimensions(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds, bool mixed = false, unsigned int currentIterNum = 0);

            unsigned int useMixedOptimizers(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, vector <unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, string* textToPrint, std::map<uint, uint> &baseNumberUpperBounds);
            void optimizeBaseNum(SingleProcessPhyloLikelihood* tl, size_t index, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, double lowerBound, double upperBound, const string &paramName, ParameterList& params, uint model, std::map<uint, uint> &baseNumberUpperBound);

            // // function working on the likelihoods vector object
            void clearVectorOfLikelihoods(size_t new_size);
            void clearVectorOfLikelihoods(size_t new_size, std::vector<SingleProcessPhyloLikelihood*> &likelihoodsVec);
            void deleteLikObject(SingleProcessPhyloLikelihood* lik_to_del);
            

            // // helper functions for optimization
            void checkLegalUseOfGradientOptimization();
            vector <string> getNonFixedParams(SingleProcessPhyloLikelihood* tl, ParameterList &allParams, map<uint, vector<int>>* fixedParams) const;
            void getBaseNumCandidates(vector <unsigned int> &baseNumCandidates, std::map<uint, uint> &baseNumUpperBounds) const;
            void fillVectorOfBaseNumCandidates(vector <unsigned int> &baseNumCandidates, unsigned int lowerBound, unsigned int upperBound) const;
            uint getMaxBaseNumAmongModels(std::map<uint, uint> baseNumberUpperBound) const;
            void getAllPossibleChrRanges(vector <unsigned int> &baseNumCandidates) const;
            //string findParameterNameInModel(string fullParameterName) const;
            //void setNewBounds(const ParameterList params, Parameter &param, map<string, pair<string, bool>> &paramPairsMap, double* lowerBound, const ChromosomeSubstitutionModel* model);

            //print functions
            void printLikParameters(SingleProcessPhyloLikelihood* lik, unsigned int optimized, string* textToPrint, const string path = "none") const;
           
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
            void updateSharedParameters(std::map<int, vector<std::pair<uint, int>>> &sharedParams, uint prevShift, uint numOfShifts) const;

            
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

    };
}
#endif // CHROMEVOL_CHROMOSOMENUMBEROPTIMIZER_H