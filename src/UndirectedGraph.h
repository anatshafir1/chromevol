#ifndef _UNDIRECTEDGRAPH_H_
#define _UNDIRECTEDGRAPH_H_

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <algorithm>

using namespace std;
namespace bpp{
    class Vertex{
        protected:
            uint nodeId_;
            vector<uint> neighbors_;
            uint root_;
        public:
            Vertex(uint nodeId):
                nodeId_(nodeId),
                neighbors_(),
                root_(nodeId)
            {}
            Vertex(uint nodeId, vector<uint> &neighbors):
                nodeId_(nodeId),
                neighbors_(neighbors),
                root_(nodeId)
            {}
            Vertex(const Vertex& v):
                nodeId_(v.nodeId_),
                neighbors_ (v.neighbors_),
                root_(v.root_)
            {}
            Vertex& operator=(const Vertex& v){
                nodeId_ = v.nodeId_;
                neighbors_ = v.neighbors_;
                root_ = v.root_;
                return *this;
            }
            Vertex* clone() const { return new Vertex(*this); }
            virtual ~Vertex(){}
        public:

            uint getId() const{
                return this->nodeId_;
            }
            const vector<uint> getNeighbors() const{
                return neighbors_;
            }
            const uint getRoot() const{
                return root_;
            }

        //protected:
            void setNeighbor(uint node){
                neighbors_.push_back(node);
            }
            void setRoot(uint root){
                root_ = root;
            }




    };

    class UndirectedGraph{
        private:
            vector<Vertex*> vertices_;
            map <pair<uint, uint>, double> edges_;
            std::map<uint, Vertex*> mapNodeIdsNodes_;
            vector<Vertex*> roots_;


        public:
            UndirectedGraph():
                vertices_(),
                edges_(),
                mapNodeIdsNodes_(),
                roots_() {}

            UndirectedGraph(Vertex* node):
                vertices_(),
                edges_(),
                mapNodeIdsNodes_(),
                roots_()
            {
                vertices_.push_back(node);
                mapNodeIdsNodes_[node->getId()] = node;
                roots_.push_back(node);
            }
            UndirectedGraph(const UndirectedGraph& graph):
                vertices_(graph.vertices_),
                edges_(graph.edges_),
                mapNodeIdsNodes_(graph.mapNodeIdsNodes_),
                roots_(graph.roots_)
            {}
            UndirectedGraph& operator=(const UndirectedGraph& graph){
                vertices_ = graph.vertices_;
                edges_ = graph.edges_;
                mapNodeIdsNodes_ = graph.mapNodeIdsNodes_;
                roots_ = graph.roots_;
                return *this;
            }
            UndirectedGraph* clone() const { return new UndirectedGraph(*this); }
            virtual ~UndirectedGraph(){
                while(vertices_.size() > 0){
                    auto v = vertices_.back(); 
                    vertices_.pop_back();
                    delete v;
     
                }
            }
            
        public:
            Vertex* getNode(uint id) const{
                return mapNodeIdsNodes_.at(id);
            }
            double getEdgeValue(pair<uint, uint> edge) const{
                return edges_.at(edge);
            }
            void addNewOrphanVertex(Vertex* node){
                vertices_.push_back(node);
                mapNodeIdsNodes_[node->getId()] = node;
                roots_.push_back(node);


            }
            std::vector<uint> getClusterRoots() const{
                std::vector<uint> clusterRoots;
                for (size_t i = 0; i < roots_.size(); i++){
                    clusterRoots.push_back(roots_[i]->getId());
                }
                return clusterRoots;
            }
            Vertex* getRoot(Vertex* v) const{
                return mapNodeIdsNodes_.at(v->getRoot());
            }
            void addEdgeBetweenTwoNodes(Vertex* node1, Vertex* node2, double edgeVal);
  
            void addEdgeBetweenTwoNodes(uint node1, uint node2, double edgeVal);

            void countNumOfVerticesInClusters(std::map<uint, vector<Vertex*>> &clusters, std::map<uint, vector<pair<uint, uint>>> &bestEdges);
            std::map<uint, bool> findFullyConnectedClusters(std::map<uint, vector<Vertex*>> &clusters);

            void printGraph() const;    
            bool hasEdges() const{
                return (vertices_.size() != roots_.size());
            }
            std::vector<uint> getNodesIds(vector<Vertex*> nodes) const{
                std::vector<uint> nodesIds;
                for (size_t i = 0; i < nodes.size(); i++){
                    nodesIds.push_back(nodes[i]->getId());
                }
                return nodesIds;
            }
        protected:
            void countNumOfVerticesInClusterRec(std::map<uint, bool> &discoveredNodes, vector<Vertex*> &cluster, Vertex* visitedNode, vector<pair<uint, uint>>* maxEdge);
            void updateRoot(Vertex* node, uint newRootValue);




    };
}
#endif // _UNDIRECTEDGRAPH_H_
