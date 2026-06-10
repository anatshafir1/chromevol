#ifndef CHROMEVOL_JOINTCONTINUOUSTRAITCHRSIMMNG_H
#define CHROMEVOL_JOINTCONTINUOUSTRAITCHRSIMMNG_H

#include "ChromEvolOptions.h"
#include "ChromosomeSubstitutionModel.h"
#include "JointEvolutionSimulator.h"
#include "LikelihoodUtils.h"
#include "ContinuousTraitModel.h"
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>

#include <string>
#include <vector>
#include <iostream>
#include <sys/stat.h>


namespace bpp{



    class JointContinuousTraitChrSimMng{
        protected:
        const PhyloTree* tree_;
        std::shared_ptr<JointEvolutionSimulator> simulator_;
        std::vector<vector<double>> gain_; // gain_[i] is the vector of gain values for the different chromosome numbers in model i
        std::vector<vector<double>> loss_; // loss_[i] is the vector of loss values for the different chromosome numbers in model i
        std::vector<vector<double>> dupl_; // dupl_[i] is the vector of duplication values for the different chromosome numbers in model i
        std::vector<vector<double>> demi_; // demi_[i] is the vector of demipolyploidy values for the different chromosome numbers in model i
        std::vector<vector<double>> baseNumR_; // baseNumR_[i] is the vector of base number rate values
        std::vector<uint> baseNum_; // baseNum_[i] is the base number for model i
        uint minChrNum_; // minChrNum_[i] is the minimum chromosome number for model i
        uint maxChrNum_; // maxChrNum_[i] is the maximum chromosome number for model i
        uint maxBaseNumTransition_; // maxBaseNumTransition_[i] is the maximum base number transition for model i
        std::vector<double> mu_; // mu_[i] is the value of the mu parameter for model i
        std::vector<vector<double>> sigma_; // sigma_[i] is the value of the sigma
        size_t numOfSimulations_;
        vector<vector<double>> simulationsResults_;
        vector<string> speciesNames_;
        vector<uint> rootChromState_;
        
        public:
            JointContinuousTraitChrSimMng(PhyloTree* tree, uint minChrNum, uint maxChrNum, uint maxBaseNumTransition): 
                tree_(tree),
                simulator_(),
                gain_(), loss_(), dupl_(), demi_(), baseNumR_(), baseNum_(), minChrNum_(minChrNum), maxChrNum_(maxChrNum), maxBaseNumTransition_(maxBaseNumTransition),
                 mu_(), sigma_(), numOfSimulations_(0), simulationsResults_(), speciesNames_(tree->getAllLeavesNames()), rootChromState_()
            {
                init();
                simulationsResults_.resize(numOfSimulations_);


            }
            JointContinuousTraitChrSimMng(const JointContinuousTraitChrSimMng &bm):
            tree_(bm.tree_),
            simulator_(bm.simulator_),
            gain_(bm.gain_),
            loss_(bm.loss_),
            dupl_(bm.dupl_),
            demi_(bm.demi_),
            baseNumR_(bm.baseNumR_),
            baseNum_(bm.baseNum_),
            minChrNum_(bm.minChrNum_),
            maxChrNum_(bm.maxChrNum_),
            maxBaseNumTransition_(bm.maxBaseNumTransition_),
            mu_(bm.mu_),
            sigma_(bm.sigma_),
            numOfSimulations_(bm.numOfSimulations_),
            simulationsResults_(bm.simulationsResults_),
            speciesNames_(bm.speciesNames_),
            rootChromState_(bm.rootChromState_)
        {}



        ~JointContinuousTraitChrSimMng() {
        }

        JointContinuousTraitChrSimMng* clone() const
        {
            return new JointContinuousTraitChrSimMng(*this);
        }

        void init();
        void simulateEvolutionForAllModels(const string &modelName);
        std::shared_ptr<ContinuousTraitModel> createTraitModel(const string& modelName, double mu, vector<double> sigma, uint minChrNum, uint maxChrNum, ChromosomeNumberDependencyFunction::FunctionType sigmaFuncType);
        void writeSimulationsResultsToFile(const string& outputFilePath);

        protected:
        void updateModelParameters(std::unordered_map<string, vector<double>>& modelParameters, const string& paramName, vector<vector<double>>& targetVector, size_t &numOfSimulations);
        void updateModelParameters(std::unordered_map<string, vector<double>>& modelParameters, const string& paramName, vector<double>& targetVector, size_t &numOfSimulations);
        void updateModelParameters(std::unordered_map<string, vector<double>>& modelParameters, const string& paramName, vector<uint>& targetVector, size_t &numOfSimulations);

        std::vector<string> getMatchingParamNames(const std::unordered_map<string, vector<double>>& modelParameters, const string& paramName);




};
}
#endif // CHROMEVOL_JOINTCONTINUOUSTRAITCHRSIMMNG_H