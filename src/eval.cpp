#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>
#include <functional>
#include <unordered_set>

#include <parallel-hashmap/phmap.h>

#include "log.h"
#include "global.h"
#include "uid-generator.h"
#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "sak.h"
#include "stream-obj.h"
#include "input-agp.h"
#include "input-filters.h"
#include "input-gfa.h"


#include "fibonacci-heap.h"
#include "alignments.h"
#include "input.h"
#include "eval.h"

void evalGFA(InSequences& InSequences, InAlignments& InAlignments) {
    
    InAlignments.buildEdgeGraph(InSequences.getHash1(), InSequences.getHash2(), InSequences.getuId());
    std::vector<std::vector<Edge>> adjEdgeList = InAlignments.getEdgeGraph();
    std::vector<InEdge>* edges = InSequences.getEdges();
    phmap::flat_hash_map<unsigned int, std::string>* idsToHeaders = InSequences.getHash2();
    std::string weight;
    
    for(InEdge& edge : *edges) {
        
        lg.verbose("Assessing edge: " + (*idsToHeaders)[edge.getsId1()] + "(" + std::to_string(edge.getsId1()) + ") " + edge.getsId1Or() + " " + (*idsToHeaders)[edge.getsId2()] + "(" + std::to_string(edge.getsId2()) + ") " + edge.getsId2Or());
        
        Edge fwEdge {edge.getsId1Or(), edge.getsId2(), edge.getsId2Or()};
        
        auto it = find(adjEdgeList.at(edge.getsId1()).begin(), adjEdgeList.at(edge.getsId1()).end(), fwEdge);
        
        if (it != adjEdgeList.at(edge.getsId1()).end()) {
            lg.verbose("Edge implied by read alignment (weight: " + std::to_string(it->weight) + ")");
            weight = std::to_string(it->weight);
        }else{
            lg.verbose("Edge not supported by read alignment (weight: 0)");
            weight = '0';
        }
            
        Tag tag {'i', "RC", weight};
        edge.appendTag(tag);
    }
}

struct Step {
    uint32_t id;
    char orientation;
};

struct Path {
    std::vector<Step> path;
    
    Path() {}
    
    Path(std::vector<Step> path) : path(path) {}
    
    void push_back(uint32_t id, char orientation) {
        path.push_back(Step{id, orientation});
    }
    
    std::unordered_set<uint32_t> pathToSet() {
        
        std::unordered_set<uint32_t> uIdSet;
        
        for (Step step : path)
            uIdSet.insert(step.id);
        return uIdSet;
    }
};

void dijkstra(InSequences& inSequences, std::vector<std::string> nodeList, std::string source, std::string destination, uint32_t maxSteps) {
    
    uint32_t steps = 0, pId = 0; // true if we reached a node in the original graph
    std::vector<uint64_t> destinations;
    FibonacciHeap<std::pair<const uint32_t,Path>*> Q; // node priority queue Q
    phmap::flat_hash_map<uint32_t,uint32_t> dist; // distance table
    phmap::flat_hash_map<std::string,uint32_t> nodes;
    inSequences.buildEdgeGraph();
    std::vector<std::vector<Edge>> &adjEdgeList = inSequences.getAdjEdgeList();
     // get the headers to uIds table to look for the header
    
    phmap::flat_hash_map<std::string, unsigned int> &headersToIds = *inSequences.getHash1();
    
    // map node names to internal ids
    nodeList.push_back(source);
    nodeList.push_back(destination);
    for (std::string node : nodeList) {
        phmap::flat_hash_map<std::string, unsigned int>::const_iterator got = headersToIds.find(node);
        if (got != headersToIds.end()) {
            nodes[node] = got->second;
        }else{
            fprintf(stderr, "Error: node not in graph (pIUd: %s)\n", node.c_str());
            exit(EXIT_FAILURE);
        }
    }
    lg.verbose("Nodelist loaded");
    dist[pId] = 0;
    Path firstPath;
    firstPath.push_back(nodes[source],'+');
    std::pair<const uint32_t,Path> *u = new std::pair<const uint32_t,Path>(pId++, firstPath);
    Q.insert(u, 0); // append source node
    lg.verbose("Starting search");
    while (Q.size() > 0 && steps < maxSteps) { // the main loop
        std::pair<const uint32_t,Path> u = *Q.extractMin(); // remove and return best segment
        InSegment &segment = inSequences.findSegmentBySUId(u.second.path.back().id);
        std::cout<<"we are at segment: "<<segment.getSeqHeader()<<std::endl;
        if (segment.getSeqHeader() == destination) {
            bool hamiltonian = true;
            std::cout<<"destination found."<<std::endl;
            std::unordered_set<uint32_t> pathNodes = u.second.pathToSet();
            for (auto n : u.second.path) {
                InSegment &segment = inSequences.findSegmentBySUId(n.id);
                std::cout<<segment.getSeqHeader()<<",";
            }
            std::cout<<std::endl;
                
            for (auto& it: nodes) {
                auto found = pathNodes.find(it.second);
                if (found == pathNodes.end()) {
                    std::cout<<"This is not a Hamiltonian path."<<std::endl;
                    dist[pId] = std::numeric_limits<uint32_t>::max();
                    hamiltonian = false;
                    break;
                }
            }
            if (hamiltonian) {
                std::cout<<"Hamiltonian path found."<<std::endl;
                break;
            }else {
                continue;
            }
        }
        uint32_t alt = dist[pId] + 1;
        for(auto v : adjEdgeList.at(segment.getuId())) {
            
            if (u.second.path.back().orientation != v.orientation0)
                continue;
            
            InSegment &nextSegment = inSequences.findSegmentBySUId(v.id);
            if (nodes.find(nextSegment.getSeqHeader()) != nodes.end()) {
                std::cout<<"inspecting segment: "<<nextSegment.getSeqHeader()<<std::endl;
                Path newPath(u.second.path);
                newPath.push_back(v.id,v.orientation1);
                std::pair<const uint32_t,Path> *u = new std::pair<const uint32_t,Path>(pId++, newPath);
                Q.insert(u, alt);
                dist[pId] = alt;
            }
        }
        ++steps;
    }
    lg.verbose("Search completed");
}
