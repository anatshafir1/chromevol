#ifndef CHROMEVOL_BROWNIANMOTIONSIMULATOR_H
#define CHROMEVOL_BROWNIANMOTIONSIMULATOR_H

#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <unordered_map>
#include <random>
#include <limits>

namespace bpp
{
/**
 * @brief The BrownianMotionSimulator class: this class simulates brownian motion trait evolution goven a phylogeny, mu, and sigma.
 *
 */ 

class BrownianMotionSimulator
{
    protected:
        double mu_;
        double sigma_;
        std::shared_ptr<PhyloTree> tree_;
        std::unordered_map<uint, double> traitStatesBranches_;
        std::unordered_map<uint, double> traitStatesInNodes_;
        int seed_;
        double minStateInLeaves_;
        double maxStateInLeaves_;
        double minStateAbsolute_;
        double maxStateAbsolute_;


    public:
        BrownianMotionSimulator(std::shared_ptr<PhyloTree> tree, double &mu, double & sigma, int seed=1):
            mu_(mu),
            sigma_(sigma),
            tree_(tree),
            traitStatesBranches_(),
            traitStatesInNodes_(),
            seed_(seed),
            minStateInLeaves_(),
            maxStateInLeaves_(),
            minStateAbsolute_(),
            maxStateAbsolute_()
        {}


        BrownianMotionSimulator(const BrownianMotionSimulator &bm):
            mu_(bm.mu_),
            sigma_(bm.sigma_),
            tree_(bm.tree_),
            traitStatesBranches_(bm.traitStatesBranches_),
            traitStatesInNodes_(bm.traitStatesInNodes_),
            seed_(bm.seed_),
            minStateInLeaves_(bm.minStateInLeaves_),
            maxStateInLeaves_(bm.maxStateInLeaves_),
            minStateAbsolute_(bm.minStateAbsolute_),
            maxStateAbsolute_(bm.maxStateAbsolute_)

        {}

        ~BrownianMotionSimulator() {}

        BrownianMotionSimulator* clone() const
        {
            return new BrownianMotionSimulator(*this);
        }


    protected:
        void simulate_(uint &nodeId, std::default_random_engine &generator);


    public:
        void simulate();
        const std::unordered_map<uint, double> getStatesInNodes() const{return traitStatesInNodes_;}
        const std::unordered_map<uint, double> getStatesInBranches() const{return traitStatesBranches_;}
        const std::shared_ptr<PhyloTree> getTree() const {return tree_;}
        const double getMinStateInData() const {return minStateInLeaves_;}
        const double getMaxStateInData() const {return maxStateInLeaves_;}
        const double getMinStateAbsolute() const {return minStateAbsolute_;}
        const double getMaxStateAbsolute() const {return maxStateAbsolute_;}
        std::unordered_map<string, double> getLeavesStates() const;


};
}
#endif // CHROMEVOL_BROWNIANMOTIONSIMULATOR_H