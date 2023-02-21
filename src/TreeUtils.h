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
#ifndef _TREEUTILS_H_
#define _TREEUTILS_H_
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>



//from bpp-seq


//from bpp-phyl
#include <Bpp/Phyl/Tree/PhyloTreeTools.h>
#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>


//standard libraries
#include <string>
#include <vector>
#include <iostream>
#include <sys/stat.h>


using namespace std;
namespace bpp{
  class TreeUtils{
    public:
      TreeUtils() {}
      virtual ~TreeUtils() {}

    public:
        static std::shared_ptr<PhyloNode> getMRCA(PhyloTree* tree, std::vector<shared_ptr<PhyloNode>> nodes);




  };

}
#endif // _TREEUTILS_H_