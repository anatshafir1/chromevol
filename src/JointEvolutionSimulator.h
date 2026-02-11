#ifndef CHROMEVOL_JOINTEVOLUTIONSIMULATOR_H
#define CHROMEVOL_JOINTEVOLUTIONSIMULATOR_H

#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include "ChromosomeSubstitutionModel.h"
#include "ContinuousTraitModel.h"
#include <Bpp/Phyl/Simulation/MutationProcess.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

namespace bpp
{
/**
 * @brief The BrownianMotionLikelihood class: this class manages all the likelihood computations of the brownian motion trait.
 * Reference: Fast likelihood calculations for comparative analyses (R. Freckleton 2012)
 *
 */ 

class JointEvolutionSimulator

{
        
    protected:
        std::shared_ptr<PhyloTree> tree_;
        const std::shared_ptr<ChromosomeBMSubstitutionModel> chromosomeModel_;
        const std::shared_ptr<ContinuousTraitModel> traitModel_;
        std::unordered_map<uint, vector<std::pair<double,double>>> simulatedEvolutionTrait_; // map of node id to a pair of (trait value, time)
        std::unordered_map<uint, MutationPath> simulatedEvolutionChromosome_; // map of node id to a pair of (chromosome number, time)
        double traitRootState_;
        size_t chromosomeRootState_;
        vector<double> cumSumChromosomeRootFrequencies_;
        std::unordered_map<uint, std::pair<double, int>> nodeStates_; // map of node id to a pair of (trait value, chrom number)
        int seed_;
        double deltaT_;
        bool successfulSimulation_;



    public:
        JointEvolutionSimulator(PhyloTree* tree, const std::shared_ptr<ChromosomeBMSubstitutionModel> chromosomeModel, const std::shared_ptr<ContinuousTraitModel> traitModel, vector<double> &chromRootFrequencies, int seed, double deltaT):
            tree_(std::shared_ptr<PhyloTree>(tree->clone())), chromosomeModel_(chromosomeModel), traitModel_(traitModel), simulatedEvolutionTrait_(), simulatedEvolutionChromosome_(),
            traitRootState_(traitModel->sampleRootState()), chromosomeRootState_(), cumSumChromosomeRootFrequencies_(), nodeStates_(), seed_(seed), deltaT_(deltaT), successfulSimulation_(true)

        {

            cumSumChromosomeRootFrequencies_ = VectorTools::cumSum(chromRootFrequencies);
            chromosomeRootState_ = RandomTools::pickFromCumSum(cumSumChromosomeRootFrequencies_);
            RandomTools::setSeed(static_cast<long>(seed_));
            

        }

        JointEvolutionSimulator(const JointEvolutionSimulator &bm):
            tree_(bm.tree_),
            chromosomeModel_(bm.chromosomeModel_),
            traitModel_(bm.traitModel_),
            simulatedEvolutionTrait_(bm.simulatedEvolutionTrait_),
            simulatedEvolutionChromosome_(bm.simulatedEvolutionChromosome_),
            traitRootState_(bm.traitRootState_),
            chromosomeRootState_(bm.chromosomeRootState_),
            cumSumChromosomeRootFrequencies_(bm.cumSumChromosomeRootFrequencies_),
            nodeStates_(bm.nodeStates_),
            seed_(bm.seed_),
            deltaT_(bm.deltaT_),
            successfulSimulation_(bm.successfulSimulation_)



        {}

        ~JointEvolutionSimulator() {
        }

        JointEvolutionSimulator* clone() const
        {
            return new JointEvolutionSimulator(*this);
        }


        void simulateEvolution(bool reminderStep = true);
        std::unordered_map<uint, vector<std::pair<double,double>>>& getSimulatedTraitEvolution(){return simulatedEvolutionTrait_;}
        std::unordered_map<uint, MutationPath>& getSimulatedChromosomeEvolution(){return simulatedEvolutionChromosome_;}
        //The chromosome states are already shifted according to the min chromosome number
        std::unordered_map<string, int> getSimulatedLeafChromosomeStates();
        std::unordered_map<string, double> getSimulatedLeafTraitStates();

        template <typename T>
        void writeAsFasta(const std::string& filename, const std::unordered_map<string, T>& dataMap){
            std::ofstream outFile(filename);
            if (!outFile.is_open()){
                throw Exception("JointEvolutionSimulator::writeAsFasta: Could not open file " + filename + " to write simulated chromosome data.");
            }
            for (const auto& pair: dataMap){
                outFile << ">" << pair.first << "\n" << pair.second << "\n";
            }
            outFile.close();
        }
        bool isSuccessfulSimulation(){return successfulSimulation_;}



    protected:

        void simulateEvolutionRec_(uint &nodeId, size_t initialChromState, double initialTraitState, bool reminderStep);
        bool isPolyploidyEvent(size_t prevState, size_t currentState) const;
        double updateTraitStateAfterPolyploidy(size_t prevState, size_t currentState, double traitState) const;
        void simulateMutationPath(uint nodeId, double &currentTraitState, size_t &currentChromState, double delta_time, MutationPath &mutationPathChrom, double timeAccumulated, size_t &prevState);







        




};
}
#endif // CHROMEVOL_JOINTEVOLUTIONSIMULATOR_H