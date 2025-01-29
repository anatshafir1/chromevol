#ifndef CHROMEVOL_BROWNIANMOTIONLIKELIHOOD_H
#define CHROMEVOL_BROWNIANMOTIONLIKELIHOOD_H

#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <unordered_map>

namespace bpp
{
/**
 * @brief The BrownianMotionLikelihood class: this class manages all the likelihood computations of the brownian motion trait.
 * Reference: Fast likelihood calculations for comparative analyses (R. Freckleton 2012)
 *
 */ 

class BrownianMotionLikelihood
{
    protected:
        double mu_;
        double sigma_;
        Eigen::ArrayXd V_; // contrast variances
        Eigen::ArrayXd u_; // trait value contrasts
        double likelihood_;
        std::shared_ptr<PhyloTree> tree_;
        const std::unordered_map<string, double> traitData_;
        std::unordered_map<uint, double> parsimonyStates_;
        bool calculatedContrasts_;
        std::unordered_map<uint, double> correctedBranchLengths_;


    public:
        BrownianMotionLikelihood(std::shared_ptr<PhyloTree> tree, const std::unordered_map<string, double> &traitData):
            mu_(0),
            sigma_(0),
            V_(),
            u_(),
            likelihood_(0),
            tree_(tree),
            traitData_(traitData),
            parsimonyStates_(),
            calculatedContrasts_(false),
            correctedBranchLengths_()
        {
            size_t numberOfTaxa = traitData_.size();
            V_ = Eigen::ArrayXd::Zero(numberOfTaxa);
            u_ = Eigen::ArrayXd::Zero(numberOfTaxa);
            auto leafNodes = tree_->getLeavesUnderNode(tree_->getNode(tree_->getRootIndex()));
            for (auto &leaf : leafNodes){
                auto leafName = leaf->getName();
                uint nodeId = tree_->getNodeIndex(leaf);
                parsimonyStates_[nodeId] = traitData_.at(leafName);
            }


        }


        BrownianMotionLikelihood(const BrownianMotionLikelihood &bm):
            mu_(bm.mu_),
            sigma_(bm.sigma_),
            V_(bm.V_),
            u_(bm.u_),
            likelihood_(bm.likelihood_),
            tree_(bm.tree_),
            traitData_(bm.traitData_),
            parsimonyStates_(bm.parsimonyStates_),
            calculatedContrasts_(bm.calculatedContrasts_),
            correctedBranchLengths_(bm.correctedBranchLengths_)

        {}

        ~BrownianMotionLikelihood() {}

        BrownianMotionLikelihood* clone() const
        {
            return new BrownianMotionLikelihood(*this);
        }


    protected:
        void calculateContrasts();
        void calculateSigmaMLE(){
            sigma_ = ((u_ * u_)/V_).sum();
            sigma_ /= static_cast<double>(u_.size());
        }
        void calculateMuMLE(){
            mu_ =  parsimonyStates_[tree_->getRootIndex()];
        }
        void updateContrastsVectors(size_t &index, uint nodeId);
        double getExpectedTraitValue(double &t1, double &t2, double &trait1, double &trait2);
        double calculateNewBranchLength(double branchLength, double &t1, double &t2);

    public:

        
        // this function assumes that the vectors u_, V_, sigma_, and mu_ are already computed
        // mu and sigma are arbitrary values (sampled during optimization)
        void calculateLikelihood(double mu, double sigma){
            if (!calculatedContrasts_){
                calculateContrasts();
            }
            double expr1 = V_.log().sum();
            double expr2 = (u_.square() / (sigma * V_)).sum();
            double expr3 = static_cast<double>(u_.size()) * std::log(2 * M_PI * sigma);
            double expr4 = ((mu - mu_) * (mu - mu_))/(sigma * V_(V_.size() - 1));
            likelihood_ = -0.5 * (expr1 + expr2 + expr3 +expr4);

        }
        const double getLikelihood(){return likelihood_;}
        const double getMuMLE(){
            if (!calculatedContrasts_){
                calculateContrasts();
            }
            return mu_;
        }
        const double getSigmaMLE(){
            if (!calculatedContrasts_){
                calculateContrasts();
            }
            return sigma_;
        }
        size_t getNumberOfParameters(){return 2;}




};
}
#endif // CHROMEVOL_BROWNIANMOTIONLIKELIHOOD_H