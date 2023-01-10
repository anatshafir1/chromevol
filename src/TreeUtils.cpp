#include "TreeUtils.h"


using namespace bpp;
using namespace std;


std::shared_ptr<PhyloNode> TreeUtils::getMRCA(PhyloTree* tree, std::vector<shared_ptr<PhyloNode>> nodes){
    shared_ptr<PhyloNode> mrca;
    vector<uint> nodesInIndices = tree->getNodeIndexes(nodes);
    size_t numOfFound = 0;
    size_t maxNumOfFound = 0;
    uint nodeIdOfMaxFound;
    vector<uint> nodesToBeFound;
    vector<uint> nodesToBeFoundInMax;


    for (size_t i = 0; i < nodes.size(); i++){
        numOfFound = 0;
        nodesToBeFound = nodesInIndices;
        if (!(tree->isLeaf(nodes[i]))){
            auto nodesUnderSubtree = tree->getSubtreeNodes(nodes[i]);
            auto subtreeNodesIndices = tree->getNodeIndexes(nodesUnderSubtree);
            for (size_t j = 0; j < subtreeNodesIndices.size(); j++){
                auto it = std::find(nodesInIndices.begin(), nodesInIndices.end(), subtreeNodesIndices[j]);

                if (it != nodesInIndices.end()){
                    nodesToBeFound.erase(std::remove(nodesToBeFound.begin(), nodesToBeFound.end(), subtreeNodesIndices[j]), nodesToBeFound.end());
                    numOfFound ++;
                    if (numOfFound == nodesInIndices.size()){
                        mrca = nodes[i];
                        return mrca;
                    }
                }
            }

        }else{
            nodesToBeFound.erase(std::remove(nodesToBeFound.begin(), nodesToBeFound.end(), nodesInIndices[i]), nodesToBeFound.end());
            numOfFound = 1;
        }

        if (numOfFound >= maxNumOfFound){
          nodeIdOfMaxFound = nodesInIndices[i];
          maxNumOfFound = numOfFound;
          nodesToBeFoundInMax = nodesToBeFound;
        }
    }
    // MRCA was not among the nodes
    // choose the one which contained the most nodes in its subtree (or if none had, one is just chosen randomely)

    uint nodeId = nodeIdOfMaxFound;
    while (maxNumOfFound < nodesInIndices.size()){
        if (nodeId == tree->getRootIndex()){
            mrca = tree->getRoot();
            break;
        }
        auto edgeIndex =  tree->getIncomingEdges(nodeId)[0]; 
        auto fatherIndex = tree->getFatherOfEdge(edgeIndex);
        auto fatherNode = tree->getNode(fatherIndex);
        auto sons = tree->getSons(fatherNode);
        for (size_t n = 0; n < sons.size(); n++){
          uint sonId = tree->getNodeIndex(sons[n]);
          if (sonId == nodeId){
            continue;
          }
          auto nodesOfSubtree = tree->getSubtreeNodes(sons[n]);
          for (size_t i = 0; i < nodesOfSubtree.size(); i++){
              auto subtreeNodeId = tree->getNodeIndex(nodesOfSubtree[i]);
              auto it = std::find(nodesToBeFoundInMax.begin(), nodesToBeFoundInMax.end(), subtreeNodeId);
              if (it != nodesToBeFoundInMax.end()){
                  nodesToBeFoundInMax.erase(std::remove(nodesToBeFoundInMax.begin(), nodesToBeFoundInMax.end(), subtreeNodeId), nodesToBeFoundInMax.end());
                  maxNumOfFound ++;
                  mrca = tree->getNode(fatherIndex);
              }
          }
        }

        nodeId = fatherIndex;

    }
    return mrca;
}
