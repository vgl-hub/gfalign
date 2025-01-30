#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>
#include <functional>

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

void dijkstra(InSequences& inSequences, std::vector<std::string> nodeList, std::string source, std::string destination, uint32_t maxSteps) {
    
    uint32_t steps = 0; // true if we reached a node in the original graph
    std::vector<uint64_t> destinations;
    FibonacciHeap<std::pair<const uint32_t,InSegment&>*> Q; // node priority queue Q
    phmap::flat_hash_map<uint32_t,uint32_t> dist; // distance table
    phmap::flat_hash_map<uint32_t,std::string> prev; // previous node
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
    dist[nodes[source]] = 0;
    InSegment &segment = inSequences.findSegmentBySUId(nodes[source]);
    std::pair<const uint32_t,InSegment&> u = std::make_pair(headersToIds[source], std::ref(segment));
    Q.insert(&u, 0); // append source node

    while (Q.size() > 0 && steps < maxSteps) { // the main loop
        std::pair<const uint32_t,InSegment&> u = *Q.extractMin(); // remove and return best vertex
            
        for(auto v : adjEdgeList.at(u.first)) {
            InSegment &segment = inSequences.findSegmentBySUId(v.id);
            if (nodes.find(segment.getSeqHeader()) != nodes.end()) {
                dist[v.id] = std::numeric_limits<uint8_t>::max();
                std::pair<const uint32_t,InSegment&> u = std::make_pair(v.id, std::ref(segment));
                Q.insert(&u, std::numeric_limits<uint8_t>::max());
                uint32_t alt = dist[u.first] + 1;
                if (alt < dist[v.id]) {
                    prev[v.id] = u.first;
                    dist[v.id] = alt;
                    Q.decreaseKey(&u, alt);
                }
            }
        }
        ++steps;
    }
//    std::deque<std::string> S;
//    std::string u = destination;
//    if (prev.find(u) != prev.end() || u == source) {
//        while (true) {
//            S.push_front(u);
//            u = prev[u];
//        }
//    }
//    for (const auto& node : S)
//        std::cout << node << " ";
//    std::cout << std::endl;
}
