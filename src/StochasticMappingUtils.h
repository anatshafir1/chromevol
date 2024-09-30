//
// File: ChromosomeNumberMng.h
// Created by: Anat Shafir
// Created on: Mon September 11 14:57 2020
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
#ifndef CHROMEVOL_STOCHASTICMAPPINGUTILS_H
#define CHROMEVOL_STOCHASTICMAPPINGUTILS_H
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>



//from bpp-seq

#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>

//from bpp-phyl
#include <Bpp/Phyl/Tree/PhyloTreeTools.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include "ChromosomeSubstitutionModel.h"


//standard libraries
#include <string>
#include <vector>
#include <iostream>
#include <sys/stat.h>


#define THRESHOLD_EXP 0.5
#define THRESHOLD_HEURISTIC 0.95
#define MAX_ITER_HEURISTICS 10
#define MAX_SIM_HEURISTICS 1000
#define BRANCH_MULTIPLIER_FACTOR 2

using namespace std;
namespace bpp{
  class StochasticMappingUtils{
    public:
      StochasticMappingUtils() {}
      virtual ~StochasticMappingUtils() {}

    public:
      static void printResultsForEachMapping(PhyloTree* tree, std::map<uint, std::map<int, double>> &expectationsPerTypeRootToLeaf, const NonHomogeneousSubstitutionProcess* NonHomoProcess, std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> &rootToLeafTransitions, std::map<uint, std::map<size_t, bool>> &presentMapping, const string &outStMappingRootToLeafPath, size_t mappingIndex, bool demiOnlyForEven);
      static void printStochasticMappingEvolutionaryPath(ChromosomeAlphabet* alphabet, std::shared_ptr<PhyloTree> stmTree, const std::map<uint, std::vector<MutationPath>> &mappings, const std::map<uint, std::vector<size_t>> &ancestralStates, size_t mappingIndex, const string &outPathPerMapping);
      static void printMappingStretchedBranches(const string &stretchedBranchesPath, std::map<uint, std::map<size_t, std::pair<double, double>>> &originalBrLenVsStretched, const std::shared_ptr<PhyloTree> tree);
      static void fixFailedMappings(PhyloTree* tree, StochasticMapping* stm, std::map<uint, std::map<size_t, std::pair<double, double>>> &originalBrLenVsStretched, size_t numOfFixingMappingIterations);
      static void printRootToLeaf(PhyloTree* tree, ChromosomeAlphabet* alphabet, std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> &rootToLeafOccurrences, std::map<uint, std::map<size_t, bool>> &presentMapping, size_t numOfMappings, const NonHomogeneousSubstitutionProcess* NonHomoProcess, const string &resultsDir, bool demiOnlyForEven);
      static void writeNanInTable(ofstream &stream);
      static void writeZeroInTable(ofstream &stream);
      // get model index for each node (i.e., the model of the father, because this is the model which applies on the son branch)
      static std::map<uint, size_t> getModelForEachBranch(PhyloTree &tree, const NonHomogeneousSubstitutionProcess &models);
      static void getModeForSons(PhyloTree &tree, uint fatherId, uint sonId, const NonHomogeneousSubstitutionProcess* NonHomoModel,  std::map<uint, size_t> &modelPerNode);
      static std::string getTypeOfTransitionStr(int transitionType);
      static std::map<int, double> getTypeForEachTransitionPerNode(std::shared_ptr<const ChromosomeSubstitutionModel> chrModel, std::map<pair<size_t, size_t>, double> &transitionsPerNode, uint nodeId, bool demiOnlyForEven);
      static bool getProbabilitiesPerType(vector<double> &probabilities, int startState, int endState, std::shared_ptr<const ChromosomeSubstitutionModel> model, bool demiOnlyForEven);
      static vector <uint> getVectorOfMapKeys(std::map<uint, vector<size_t>> &mapOfVectors);


  };

}
#endif // CHROMEVOL_STOCHASTICMAPPINGUTILS_H