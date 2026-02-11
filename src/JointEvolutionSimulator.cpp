#include "JointEvolutionSimulator.h"

using namespace std;
using namespace bpp;

void JointEvolutionSimulator::simulateEvolution(bool reminderStep){
    uint rootIndex = tree_->getRootIndex();
    simulateEvolutionRec_(rootIndex, chromosomeRootState_, traitRootState_, reminderStep);

}

void JointEvolutionSimulator::simulateEvolutionRec_(uint &nodeId, size_t initialChromState, double initialTraitState, bool reminderStep){
    if (tree_->isLeaf(nodeId)){
        return;
    }
    auto sons = tree_->getSons(nodeId);
    for (auto &sonId: sons){
        double currentTraitState = initialTraitState;
        size_t currentChromState = initialChromState;
        size_t prevState = initialChromState;
        auto branch = tree_->getEdgeToFather(sonId);
        double branchLength = branch->getLength();
        MutationPath mutationPathChrom = MutationPath(chromosomeModel_->getAlphabet(), initialChromState, branchLength);
        size_t numSegments = static_cast<size_t>(branchLength / deltaT_);
        double remainingTime = branchLength - (static_cast<double>(numSegments) * deltaT_);
        double delta_time = deltaT_;
        if (remainingTime > 0.0){
            numSegments += 1;
        }

        for (size_t seg = 0; seg < numSegments; seg++){
            if ((seg == numSegments-1) && (remainingTime > 0.0)){
                delta_time = remainingTime;          
            }
            double timeAccumulated = 0.0;
            simulateMutationPath(sonId, currentTraitState, currentChromState, delta_time, mutationPathChrom, timeAccumulated, prevState);

        }
        nodeStates_[sonId] = std::make_pair(currentTraitState, static_cast<int>(currentChromState) + chromosomeModel_->getMin());

        
        simulateEvolutionRec_(sonId, currentChromState, currentTraitState, reminderStep);
        simulatedEvolutionChromosome_.emplace(sonId, mutationPathChrom);
        
    }


}
void JointEvolutionSimulator::simulateMutationPath(uint nodeId, double &currentTraitState, size_t &currentChromState, double delta_time, MutationPath &mutationPathChrom, double timeAccumulated, size_t &prevState){
    if (timeAccumulated >= delta_time){
        return;
    }
    auto chromModel = std::shared_ptr<ChromosomeBMSubstitutionModel>(chromosomeModel_->getModelAccordingToState(chromosomeModel_, currentTraitState));
    auto simulator = SimpleMutationProcess(chromModel.get());
    MutationPath mutationPath = simulator.detailedEvolve(currentChromState, delta_time-timeAccumulated);
    vector<double> times = mutationPath.getTimes();
    vector<size_t> states = mutationPath.getStates();
    bool polyploidy = false;

            
    for (size_t i = 0; i < times.size(); i++){
        if (states[i] == chromosomeModel_->getMax() - chromosomeModel_->getMin()){
            successfulSimulation_ = false;
            return;
        }
        mutationPathChrom.addEvent(states[i], times[i]);
        currentChromState = states[i];
        timeAccumulated += times[i];
        polyploidy = isPolyploidyEvent(prevState, states[i]);
        currentTraitState = traitModel_->simulate(times[i], currentTraitState, static_cast<double>(prevState) + chromosomeModel_->getMin());
        simulatedEvolutionTrait_[nodeId].push_back(std::make_pair(currentTraitState, times[i]));
        if (polyploidy){
            currentTraitState = updateTraitStateAfterPolyploidy(prevState, states[i], currentTraitState);
        }
        prevState = states[i];
        if (polyploidy){
            break;
        }
    }
    if (polyploidy){
        simulateMutationPath(nodeId, currentTraitState, currentChromState, delta_time, mutationPathChrom, timeAccumulated, prevState);
    }
    if (timeAccumulated >= delta_time){
        return;
    }

    currentTraitState = traitModel_->simulate(delta_time-timeAccumulated, currentTraitState, static_cast<double>(currentChromState) + chromosomeModel_->getMin());
    simulatedEvolutionTrait_[nodeId].push_back(std::make_pair(currentTraitState, delta_time-timeAccumulated));


}


// what if # segments which is branch_length / deltaT is not an integer?
// then we can do floor and then add one more step with the remaining time
// or we can just do ceil and reduce deltaT accordingly, so that the number of segments multiplied by the new deltaT will give the branch length

bool JointEvolutionSimulator::isPolyploidyEvent(size_t prevState, size_t currentState) const{
    
    if (currentState == prevState  + 1){
        return false;
    }else if (currentState == prevState -1){
        return false;
    }else if (currentState == prevState){
        return false;
    }
    return true;
}

double JointEvolutionSimulator::updateTraitStateAfterPolyploidy(size_t prevState, size_t currentState, double traitState) const{
    int minChrNum = chromosomeModel_->getMin();
    int chromStatePrev = static_cast<int>(prevState) + minChrNum;
    int chromStateCurr = static_cast<int>(currentState) + minChrNum;
    double newTraitState = traitState * (static_cast<double>(chromStateCurr) / static_cast<double>(chromStatePrev));
    return newTraitState;

}



std::unordered_map<string, int> JointEvolutionSimulator::getSimulatedLeafChromosomeStates(){
    std::unordered_map<string, int> leafChromStates;
    auto treeLeaves = tree_->getAllLeaves();
    for (auto &leaf: treeLeaves){
        uint leafId = tree_->getNodeIndex(leaf);
        auto leafName = leaf->getName();
        leafChromStates[leafName] = nodeStates_[leafId].second;
    }
    return leafChromStates;
}
std::unordered_map<string, double> JointEvolutionSimulator::getSimulatedLeafTraitStates(){
    std::unordered_map<string, double> leafTraitStates;
    auto treeLeaves = tree_->getAllLeaves();
    for (auto &leaf: treeLeaves){
        uint leafId = tree_->getNodeIndex(leaf);
        auto leafName = leaf->getName();
        leafTraitStates[leafName] = nodeStates_[leafId].first;
    }
    return leafTraitStates;
}

