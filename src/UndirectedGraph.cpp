#include "UndirectedGraph.h"


using namespace bpp;
using namespace std;

void UndirectedGraph::printGraph() const{
    for (size_t i = 0; i < vertices_.size(); i++){
        auto node = vertices_[i];
        std::cout << "#Node is: "<< node->getId() << std::endl;
        auto neighbors  = node->getNeighbors();
        for (size_t j = 0; j < neighbors.size(); j++){
            pair <uint, uint> edge(node->getId(), neighbors[j]);
            std::cout << "\t" << "edge: " << node->getId() << "->" << neighbors[j] << " is: "<< edges_.at(edge) << std::endl;
        }
    }
}
void UndirectedGraph::addEdgeBetweenTwoNodes(Vertex* node1, Vertex* node2, double edgeVal){
    node1->setNeighbor(node2->getId());
    node2->setNeighbor(node1->getId());
    std::pair<uint, uint> edge1(node2->getId(), node1->getId());
    std::pair<uint, uint> edge2(node1->getId(), node2->getId());
    edges_[edge1] = edgeVal;
    edges_[edge2] = edgeVal;
    Vertex* rootToRemove;
    Vertex* rootToUpdate;
    Vertex* nodeToUpdate;
    if (node1->getRoot() == node2->getRoot()){
        return;
    }
    if (node1->getRoot() > node2->getRoot()){
        rootToRemove = getRoot(node1);
        rootToUpdate = getRoot(node2);
        nodeToUpdate = node1;

    }else{
        rootToRemove = getRoot(node2);
        rootToUpdate = getRoot(node1);
        nodeToUpdate = node2;

    }
    updateRoot(nodeToUpdate, rootToUpdate->getId());
    roots_.erase(std::remove(roots_.begin(), roots_.end(), rootToRemove), roots_.end());


}
void UndirectedGraph::addEdgeBetweenTwoNodes(uint node1, uint node2, double edgeVal){
    addEdgeBetweenTwoNodes(getNode(node1), getNode(node2), edgeVal);
}

void UndirectedGraph::updateRoot(Vertex* node, uint newRootValue){
    node->setRoot(newRootValue);
    auto neighbors = node->getNeighbors();
    for (size_t i = 0; i < neighbors.size();i++){
        auto neighbor = getNode(neighbors[i]);
        if (getRoot(neighbor)->getId() != newRootValue){
            updateRoot(neighbor, newRootValue);

        }
        
    }
}

void UndirectedGraph::countNumOfVerticesInClusterRec(std::map<uint, bool> &discoveredNodes, vector<Vertex*> &cluster, Vertex* visitedNode, vector<pair<uint, uint>> *maxEdges){
    discoveredNodes[visitedNode->getId()] = true;
    cluster.push_back(visitedNode);
    auto neighbors = visitedNode->getNeighbors();

    for (size_t i = 0; i < neighbors.size(); i++){
        if (!(discoveredNodes[neighbors[i]])){
            countNumOfVerticesInClusterRec(discoveredNodes, cluster, getNode(neighbors[i]), maxEdges);
            if (maxEdges){
                pair<uint,uint> edge(visitedNode->getId(), neighbors[i]);
                if (maxEdges->size() > 0){
                    auto bestEdge = maxEdges->back();
                    auto bestEdgeValue = edges_[bestEdge];
                    auto candidateEdgeValue = edges_[edge];
                    if (candidateEdgeValue < bestEdgeValue){
                        while(maxEdges->size() > 0){
                            maxEdges->pop_back();
                        }
                        maxEdges->push_back(edge);
                    }else{
                        if (candidateEdgeValue == bestEdgeValue){
                            maxEdges->push_back(edge);
                        }
                    }

                }else{
                    maxEdges->push_back(edge);
                }
            }

        }
    }

    return;

}

void UndirectedGraph::countNumOfVerticesInClusters(std::map<uint, vector<Vertex*>> &clusters, std::map<uint, vector<std::pair<uint, uint>>> &maxEdges){
    std::map<uint, bool> discoveredNodes;
    
    for (size_t k = 0; k < vertices_.size(); k++){
        discoveredNodes[vertices_[k]->getId()] = false;
    }
    for (size_t i = 0; i < roots_.size(); i++){
        countNumOfVerticesInClusterRec(discoveredNodes, clusters[roots_[i]->getId()], roots_[i], &(maxEdges[roots_[i]->getId()]));

    }
    
    return;

}
std::map<uint, bool> UndirectedGraph::findFullyConnectedClusters(std::map<uint, vector<Vertex*>> &clusters){
    std::map<uint, bool> fullyConnectedClusters;
    for (size_t i = 0; i < roots_.size(); i++){
        auto clusterVertices = clusters[roots_[i]->getId()];
        size_t numOfVertices = clusterVertices.size();
        bool fullyConnected = true;
        for (size_t j = 0; j < clusterVertices.size(); j++){
            auto neighbors = clusterVertices[j]->getNeighbors();
            if (neighbors.size() != numOfVertices-1){
                fullyConnected = false;
                break;
            }
        }
        fullyConnectedClusters[roots_[i]->getId()] = fullyConnected;

    }
    return fullyConnectedClusters;

}