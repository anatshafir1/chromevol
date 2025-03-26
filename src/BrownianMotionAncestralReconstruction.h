#ifndef CHROMEVOL_BROWNIANMOTIONANCESTRALRECONSTRUCTION_H
#define CHROMEVOL_BROWNIANMOTIONANCESTRALRECONSTRUCTION_H

#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <unordered_map>

namespace bpp
{
/**
 * @brief The BrownianMotionAncestralReconstruction class: this class is needed to compute the maximum likelihood ancestral states given arbitratry values of mu and sigma
 *  Reference: A Linear-Time Algorithm for Gaussian and non-Gaussian Trait Evolution Models (Ho and Ane 2014)
 *  The adjustments to the arbitrary values of mu and sigma (that are not necessarily MLE) were made by myself.
 *
 */ 

class BrownianMotionAncestralReconstruction
{
    protected:
        const double muMLE_;
        const double sigmaMLE_;
        double mu_;
        double sigma_;
        std::unordered_map<string, double> leafTraitStates_; // this object can also be modified during the algorithm
        std::shared_ptr<PhyloTree> tree_; // this object should be passed cloned, because otherwise the tree will be modified!
        std::unordered_map<uint, double> ancestralTraitStates_;
        std::unordered_map<uint, double> p_;

    public:
        BrownianMotionAncestralReconstruction(double muMLE, double sigmaMLE, double mu, double sigma, std::unordered_map<string, double> leafTraitStates, std::shared_ptr<PhyloTree> tree):
            muMLE_(muMLE),
            sigmaMLE_(sigmaMLE),
            mu_(mu),
            sigma_(sigma),
            leafTraitStates_(leafTraitStates),
            tree_(tree),
            ancestralTraitStates_(),
            p_()
        {}


        BrownianMotionAncestralReconstruction(const BrownianMotionAncestralReconstruction &bm):
            muMLE_(bm.muMLE_),
            sigmaMLE_(bm.sigmaMLE_),
            mu_(bm.mu_),
            sigma_(bm.sigma_),
            leafTraitStates_(bm.leafTraitStates_),
            tree_(bm.tree_),
            ancestralTraitStates_(bm.ancestralTraitStates_),
            p_(bm.p_)

        {}

        ~BrownianMotionAncestralReconstruction() {}

        BrownianMotionAncestralReconstruction* clone() const
        {
            return new BrownianMotionAncestralReconstruction(*this);
        }

        void reconstructAncestralStates();
        const std::unordered_map<uint, double> getAncestralStates() const{return ancestralTraitStates_;}
        void claculateNodeDepths(std::unordered_map<uint, double> &nodeDepths);

    protected:

        void reconstructAncestralStatesPostOrder(uint rootId);
        void reconstructAncestralStatesPreOrder(uint rootId);
        void addLeavesToAncestralStates();
        void correctMLStatesAccordingToMu();
        void claculateNodeDepthsRec(std::unordered_map<uint, double> &nodeDepths, uint nodeId);
        //void claculateNodeDepths(std::unordered_map<uint, double> &nodeDepths);
        //void revertLeafStates();



};
}
#endif // CHROMEVOL_BROWNIANMOTIONANCESTRALRECONSTRUCTION_H
