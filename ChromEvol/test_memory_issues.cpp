#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/JointPhyloLikelihood.h>
#include <Bpp/Phyl/Model/Character/CharacterSubstitutionModel.h>
#include <Bpp/Phyl/Model/Character/SingleRateModel.h>
#include <Bpp/Phyl/Model/Character/RatePerEntryModel.h>
#include <Bpp/Phyl/Model/Character/RatePerExitModel.h>
#include <Bpp/Phyl/Model/Character/RatePerPairSymModel.h>
#include <Bpp/Phyl/Model/Character/RatePerPairModel.h>
#include "ChromosomeNumberOptimizer.h"
#include "JointTraitChromosomeLikelihood.h"
#include <Bpp/Numeric/Random/RandomTools.h>
#include <iostream>
#include <fstream>
#include <regex>

using namespace bpp;
using namespace std;

void testMemoryTree(){
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((S1:0.1,S2:0.1):0.3,S3:0.4):0.2,(S4:0.3,S5:0.3):0.3);", false, "", false, false));
  std::cout << "just created a tree" << std::endl;

  
}
int main(){
    testMemoryTree();
    return 0;

}