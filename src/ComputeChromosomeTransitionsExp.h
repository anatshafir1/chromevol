//
// File: ComputeChromosomeTransitionsExp.h
// Created by: Anat Shafir
// Created on: Mon August 23 16:48 2020
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
#ifndef _COMPUTECHROMOSOMETRANSITIONSEXP_H_
#define _COMPUTECHROMOSOMETRANSITIONSEXP_H_


#include "ChromosomeSubstitutionModel.h"
#include "Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h"
#include "Bpp/Phyl/Tree/PhyloTree.h"
#include "Bpp/Phyl/Tree/PhyloTreeTools.h"
#include "ChromEvolOptions.h"
#include "TreeUtils.h"
#include "Bpp/Phyl/Mapping/StochasticMapping.h"
#include "StochasticMappingUtils.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>

// From Seqlib:
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <vector>
#include <map>
#include <utility>
#include <string>


using namespace std;

namespace bpp
{
    class ComputeChromosomeTransitionsExp{
        // a class used to calculate the expectation of transitions per each type of 
        // transtion along the tree or per branch
        typedef pair<uint, PhyloBranch> Branch;
        private:
            // stores the joint probabilities of father and son for each branch.
            // key = node (the son of the branch). Value is a matrix represented by 
            // vector of double vectors. The first dimension refers to the son, and the second
            // to the father
            map<uint, map<size_t, VVdouble>> jointProbabilitiesFatherSon_;

            const PhyloTree* tree_;
            const NonHomogeneousSubstitutionProcess* model_;
            const ChromosomeAlphabet* alphabet_;
            // The branches on which the chromsome number changes are simulated
            vector<vector<Branch>> branchOrder_;

            // A map which stores for each node the simulated ancestral terminals, such that
            // each pair of ancestral terminals also serves as a key, where the value is a pair,
            // where the first element is the number of occurences of a given pair of ancestral terminals,
            // and the second element is a vector of expectations per each type of transition (the index represents a type).
            map <uint, map <pair<int, int>, std::pair<int, Vdouble>>> branchTransitionsExp_;

            // A map which stores for each node a map where the key is the type of chnage, and the value is the expectation of this type of change taking
            // into account all the encountered terminals. 
            map <uint, map <int, double>> expNumOfChangesPerBranch_;

            // A map where the key is a node id, and the value is the exected number of changes over all types and ancestral terminals.
            map <uint, double> expNumOfChanges_;

            int jumpTypeMethod_;    // which function to use for type classification- 0 if deterministic, 1 if probabilistic

            // A map where the key represents the transitions from state i to j, and the value is a map,
            // where the key is the type of transition, and the value is the probability that the given transition corresponds to that
            // type of transition.
            map <pair<int, int>, map<int, double>> stateJumpTypeProb_;

            /***********************************/
            // Internal functions
            /**********************************/ 

            // After the standard procedure of expextations computations, check if there any branches with not enough simulated chnages
            // i.e., the accounted ancestral terminals cover at least 95% of the possible ancestral pairs.
            // The function updates the already accounted ancestral terminals, so that they will be not taken into 
            // consideration in the second round of simulations.
            // Returns the probability covered by the accounted for ancestral terminals.
            double isNeededHeuristics(uint nodeId, map <uint, vector<pair<int,int>>>* unAccountedNodesAndTerminals);

            // Returns the cumulative probability of the accounted changes, and updates terminalsToAccount with the accounted for ancestral terminals if provided.
            double getCumulativeProbability(uint nodeId, vector <pair<int, int>>* terminalsToAccount = 0);
            
             
            void updateNumNonAccountedBranches(map <uint, vector<pair<int,int>>>* unAccountedNodesAndTerminals, int iteration, size_t modelIndex, const string FilePath);
            void updateBranchLengths(int initState, int iteration, size_t modelIndex, map <int, double>* ratesPerState);
            void getPosteriorAndExpForNonAccountedFor(map <uint, vector<pair<int, int>>>& nonAccountedForBranchesFromFirstRun);
            void computeExpPerTypeHeuristics(map <uint, vector<pair<int, int>>>& nonAccountedForBranchesFromFirstRun);
            bool isMaxStateValid(int prevState, const ChromosomeSubstitutionModel* model) const;
            
        public:
            ComputeChromosomeTransitionsExp(const std::shared_ptr<NonHomogeneousSubstitutionProcess> model,  const PhyloTree* tree, const ChromosomeAlphabet* alphabet, map<uint, map<size_t, VVdouble>>& jointProbabilitiesFatherSon, int method = 0)
            :jointProbabilitiesFatherSon_(jointProbabilitiesFatherSon), tree_(tree), model_(model.get()), alphabet_(alphabet),
            //waitingTimes_(), jumpProbs_(), 
            branchOrder_(), 
            //ancestralTerminalsCounts_(), 
            branchTransitionsExp_(), expNumOfChangesPerBranch_(), expNumOfChanges_(), jumpTypeMethod_(method), stateJumpTypeProb_(){}

            ComputeChromosomeTransitionsExp(const ComputeChromosomeTransitionsExp& exp):
                jointProbabilitiesFatherSon_(exp.jointProbabilitiesFatherSon_),
                tree_ (exp.tree_),
                model_ (exp.model_),
                alphabet_(exp.alphabet_),
                //waitingTimes_(exp.waitingTimes_),
                //jumpProbs_(exp.jumpProbs_),
                branchOrder_(exp.branchOrder_),
                //ancestralTerminalsCounts_(exp.ancestralTerminalsCounts_),
                branchTransitionsExp_(exp.branchTransitionsExp_),
                expNumOfChangesPerBranch_ (exp.expNumOfChangesPerBranch_),
                expNumOfChanges_ (exp.expNumOfChanges_),
                jumpTypeMethod_(exp.jumpTypeMethod_),
                stateJumpTypeProb_(exp.stateJumpTypeProb_)

            {}
            ComputeChromosomeTransitionsExp& operator=(const ComputeChromosomeTransitionsExp& exp){
                jointProbabilitiesFatherSon_ = exp.jointProbabilitiesFatherSon_;
                tree_ = exp.tree_;
                model_ = exp.model_;
                alphabet_ = exp.alphabet_;
                //waitingTimes_ = exp.waitingTimes_;
                //jumpProbs_ = exp.jumpProbs_;
                branchOrder_ = exp.branchOrder_;
                //ancestralTerminalsCounts_ = exp.ancestralTerminalsCounts_;
                branchTransitionsExp_ = exp.branchTransitionsExp_;
                expNumOfChangesPerBranch_ = exp.expNumOfChangesPerBranch_;
                expNumOfChanges_ = exp.expNumOfChanges_;
                jumpTypeMethod_ = exp.jumpTypeMethod_;
                stateJumpTypeProb_ = exp.stateJumpTypeProb_;
                return *this;
            }
            ComputeChromosomeTransitionsExp* clone() const { return new ComputeChromosomeTransitionsExp(*this); }
            virtual ~ComputeChromosomeTransitionsExp(){};
            void init();
            void computeExpectationOfChangePerBranch(uint nodeId, VVdouble &jointProbFatherNode, int transitionType);
            //ChromosomeSubstitutionModel::typeOfTransition getTypeOfTransition(int startState, int endState);
            //more sophisticated function: if there is an overlap between different transition types-> the chosen state is sampled according to probabilities
            //ChromosomeSubstitutionModel::typeOfTransition getTypeOfTransitionWithProb(int startState, int endState);


            void computeExpectationPerType();
            void printResults(const string path = "none");
            PhyloTree* getResultTree();
            // // from previous used class
            void runIteration(int state, size_t modelIndex, map <uint, vector<pair<int,int>>>* unAccountedNodesAndTerminals = 0);
            void computeExpectationAndPosterior();
            void runSimulations(int numOfSimulations);
            static bool compareBranches(Branch& edge1, Branch& edge2);//sorting function to sort the branches in ascending order of length
            int getRandomState(int currentState, const ChromosomeSubstitutionModel* model);
            double getExpectation(uint nodeId, int startAncestral, int endAncestral, int typeOfChange);
            void updateMapOfJumps(int startState, int endState, const ChromosomeSubstitutionModel* model);
            void updateExpectationsPerBranch(uint nodeId, pair<int, int> ancestralTerminals, pair<int, int> jumpStates);
            void runHeuristics(const string FilePath = "none");
            

            //*** *** ***
            // Temporarily include function for dealing with chromosome number model related stochastic mapping
            // as inferred from the StochasticMapping class. these functions are needed mainly to test the method
            //*** *** ***

            // find the expectations of each transition type
            static std::map<int, double> getExpectationsPerType(const NonHomogeneousSubstitutionProcess* NonHomoProcess, PhyloTree &tree, std::map<uint, std::map<pair<size_t, size_t>, double>> &expectationsPerNode);
            
            


            
            
            
            


    };
}
#endif // _COMPUTECHROMOSOMETRANSITIONSEXP_H_