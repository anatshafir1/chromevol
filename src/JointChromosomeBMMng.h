#ifndef CHROMEVOL_JOINTCHROMOSOMEBMMNG_H
#define CHROMEVOL_JOINTCHROMOSOMEBMMNG_H

#include <Bpp/Seq/Container/VectorSiteContainer.h>
//#include "TraitOptimizer.h"
#include "ChromosomeNumberOptimizer.h"
#include "JointPhyloChromosomeBMLikelihood.h"
#include "ChromosomeSubstitutionModel.h"
#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace bpp
{
/**
 * @brief The BrownianMotionLikelihood class: this class manages all the likelihood computations of the brownian motion trait.
 * Reference: Fast likelihood calculations for comparative analyses (R. Freckleton 2012)
 *
 */ 

class JointChromosomeBMMng:
    public BaseNumberOptimizer//,
    //public TraitOptimizer
{
    protected:
    vector <JointPhyloChromosomeBMLikelihood*> vectorOfJointLikelohoods_;
    std::shared_ptr<PhyloTree> tree_;
    std::unordered_map<string, double> traitData_;
    const VectorSiteContainer* chromosomeData_;
    std::vector <unsigned int> baseNumberCandidates_;
    double tolerance_;
    uint numOfIterationsPerStep_;
    vector<unsigned int> numOfPoints_;
    vector<unsigned int> numOfIterations_;
    double minTraitValue_;
    double maxTraitValue_;
    double minTraitStateInData_;
    double maxTraitStateInData_;
    vector<string> traitParamNames_;
    vector<string> chromosomeParamNames_;

    
    protected:
    void clearVectorOfLikelihoods(vector<JointPhyloChromosomeBMLikelihood*> &vectorOLikelihoods, size_t new_size);
    void setTraitData(const string &filePath);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> initHeterogeneousModel(std::shared_ptr<ChromosomeBMSubstitutionModel> chrModel, std::shared_ptr<ParametrizablePhyloTree> tree) const;
    JointPhyloChromosomeBMLikelihood* initJointBMModel(std::shared_ptr<ParametrizablePhyloTree> parTree, std::shared_ptr<NonHomogeneousSubstitutionProcess> modelSet);
    void optimizePointInCycle(double tol, uint numOfIterations, uint numberOfIterationsPerOneOptimization, JointPhyloChromosomeBMLikelihood* likPhyloObj);
    void separateTraitAndChromosomeParamNames(JointPhyloChromosomeBMLikelihood* lik, vector<string> &traitParamNames, vector<string> &chromosomeParamNames) const;
    void optimizeChromosomeModel(double tol, uint numOfIterations, std::vector<string> &chromosomeParamNames, std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>> &paramNameAndType, JointPhyloChromosomeBMLikelihood* likPhyloObj);
    double getRoughSigmaEstimator();
    void optimizeTraitModel(JointPhyloChromosomeBMLikelihood* lik, double tol, uint numOfIterations, std::vector<string> &traitParamNames, std::vector<string> &fixedParamsTrait);
    void printLikelihoodParameters(JointPhyloChromosomeBMLikelihood* lik, unsigned int optimized, vector<string> paramsNames) const;
    void determineMinMaxStates();
    string getParamWithoutSuffixAndPrefix(string &paramName);



    public:
        JointChromosomeBMMng(PhyloTree* tree, const string &traitDataPath, std::map<uint, uint> &baseNumberUpperBound, const VectorSiteContainer* chromosomeVsc):
            BaseNumberOptimizer(chromosomeVsc, baseNumberUpperBound),
            vectorOfJointLikelohoods_(),
            tree_(std::shared_ptr<PhyloTree>(tree->clone())),
            traitData_(),
            chromosomeData_(chromosomeVsc),
            baseNumberCandidates_(),
            tolerance_(),
            numOfIterationsPerStep_(),
            numOfPoints_(),
            numOfIterations_(),
            minTraitValue_(),
            maxTraitValue_(),
            minTraitStateInData_(),
            maxTraitStateInData_(),
            traitParamNames_(),
            chromosomeParamNames_()
        {
            setTraitData(traitDataPath);
            determineMinMaxStates();
            getBaseNumCandidates(baseNumberCandidates_, baseNumberUpperBound);

        }
        void initJointPhyloChromosomeBMLikelihood(ChromosomeAlphabet* alphabet, pair<int, std::map<int, std::vector<double>>> &modelsParams, ChromosomeSubstitutionModel::rootFreqType freqType,
            vector<int> &rateChangeType,
            bool demiOnlyForEven,
            double sigma, double mu, double estimatedSigma, double parsimonyBound, bool random=false);

        void initMultipleLikelihoodPoints(ChromosomeAlphabet* alphabet, 
            pair<int, std::map<int, std::vector<double>>> &modelsParams, 
            ChromosomeSubstitutionModel::rootFreqType freqType,
            vector<int> &rateChangeType,
            bool demiOnlyForEven,
            double sigma,
            double mu,
            double parsimonyBound);


        JointChromosomeBMMng(const JointChromosomeBMMng &bm):
            BaseNumberOptimizer(bm),
            vectorOfJointLikelohoods_(bm.vectorOfJointLikelohoods_),
            tree_(bm.tree_),
            traitData_(bm.traitData_),
            chromosomeData_(bm.chromosomeData_),
            baseNumberCandidates_(bm.baseNumberCandidates_),
            tolerance_(bm.tolerance_),
            numOfIterationsPerStep_(bm.numOfIterationsPerStep_),
            numOfPoints_(bm.numOfPoints_),
            numOfIterations_(bm.numOfIterations_),
            minTraitValue_(bm.minTraitValue_),
            maxTraitValue_(bm.maxTraitValue_),
            minTraitStateInData_(bm.minTraitStateInData_),
            maxTraitStateInData_(bm.maxTraitStateInData_),
            traitParamNames_(bm.traitParamNames_),
            chromosomeParamNames_(bm.chromosomeParamNames_)


        {}

        ~JointChromosomeBMMng() {
            clearVectorOfLikelihoods(vectorOfJointLikelohoods_, 0);
        }

        JointChromosomeBMMng* clone() const
        {
            return new JointChromosomeBMMng(*this);
        }
        void initOptimizerSettings(string &baseNumOptimizationMethod, double tolerance, uint numberOfIterationsPerStep, vector<unsigned int> numOfPoints, vector<unsigned int> numOfIterations){
            setBaseNumberOptimizationMethod(baseNumOptimizationMethod);
            tolerance_ = tolerance;
            numOfIterationsPerStep_ = numberOfIterationsPerStep;
            numOfPoints_ = numOfPoints;
            numOfIterations_ = numOfIterations;

        }
        static bool compareJointLikValues(JointPhyloChromosomeBMLikelihood* lik1, JointPhyloChromosomeBMLikelihood* lik2){
            return (lik1->getValue() < lik2->getValue());
        }
        void optimizeJointLikelihood();
        const std::unordered_map<string, double> getTraitData(){return traitData_;}

        double getFinalLikelihood(){return vectorOfJointLikelohoods_[0]->getValue();}
        const JointPhyloChromosomeBMLikelihood* getBestLikelihoodObject(){return vectorOfJointLikelohoods_[0];}
        size_t getNumberOfParameters() const;
        vector<string> getModelParameters() const{
            vector<string> paramNames;
            for(auto paramName : chromosomeParamNames_){
                paramNames.push_back(paramName);

            }
            for (auto paramName : traitParamNames_){
                paramNames.push_back(paramName);
            }
            return paramNames;
        }
        std::shared_ptr<NonHomogeneousSubstitutionProcess> setSubstitutionProcess(ValueRef <Eigen::RowVectorXd> rootFreqs, SingleProcessPhyloLikelihood* chromosomeLik, const JointPhyloChromosomeBMLikelihood* lik, std::shared_ptr<ParametrizablePhyloTree> parTree = nullptr);
        std::shared_ptr<LikelihoodCalculationSingleProcess> setFixedJointPhyloBMLikelihoodForMLAncestralReconstruction();
        std::unordered_map<uint, double> getBranchStatesForEachNode();
        const std::shared_ptr<PhyloTree> getTree(){return tree_;}
        static std::shared_ptr<NonHomogeneousSubstitutionProcess> setHeterogeneousTraitDependentModel(std::shared_ptr<ChromosomeBMSubstitutionModel> chrModel, std::shared_ptr<ParametrizablePhyloTree> parTree, const std::unordered_map<uint, double> &statesInBranches, const PhyloTree* tree, std::shared_ptr<FrequencySet> &freqs);



        




};
}
#endif // CHROMEVOL_JOINTCHROMOSOMEBMMNG_H