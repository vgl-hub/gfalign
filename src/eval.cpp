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
#include "nodetable.h"
#include "alignments.h"
#include "input-gfalign.h"
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

bool evaluatePath(const Path &path, InSequences &inSequences, std::vector<Path> alignmentPaths, phmap::flat_hash_map<uint32_t,uint32_t> &dist, NodeTable &nodeTable, const uint32_t pId, const uint32_t destinationId) {
    
    Step lastStep = path.path.back();
    InSegment &segment = inSequences.findSegmentBySUId(lastStep.id);
    lg.verbose("We are at segment: " + segment.getSeqHeader() + lastStep.orientation);
    if (lastStep.id != destinationId)
        return false;
    
    lg.verbose("Destination found.");
    std::unordered_set<uint32_t> pathNodes = path.pathToSet();
    std::vector<std::string> uniques;
    for (auto n : path.path) {
        InSegment &segment = inSequences.findSegmentBySUId(n.id);
        std::cout<<segment.getSeqHeader()<<n.orientation<<",";
        uniques.push_back(segment.getSeqHeader());
    }
    std::sort(uniques.begin(), uniques.end());
    auto last = std::unique(uniques.begin(), uniques.end());
    uniques.erase(last, uniques.end());
    std::cout<<" "<<uniques.size()<<std::endl;
	
	int dp[MAX_N][MAX_N];
	for (Path alignmentPath : alignmentPaths)
		alignPaths(1, -1, -1, path, alignmentPath, dp);
    
    bool hamiltonian = true;
    for (auto& it: nodeTable.records) {
        auto found = pathNodes.find(it.second.uId);
        if (found == pathNodes.end()) {
            lg.verbose("This is not a Hamiltonian path.");
            dist[pId] = std::numeric_limits<uint32_t>::max();
            hamiltonian = false;
            break;
        }
    }
    if (hamiltonian)
        lg.verbose("Hamiltonian path found.");
    return true;
}

void dijkstra(InSequences &inSequences, InAlignments& inAlignments, std::string nodeFile, std::string source, std::string destination, uint32_t maxSteps) {
    
    uint32_t steps = 0, pId = 0; // true if we reached a node in the original graph
    std::vector<uint64_t> destinations;
    FibonacciHeap<std::pair<const uint32_t,Path>*> Q; // node priority queue Q
    phmap::flat_hash_map<uint32_t,uint32_t> dist; // distance table
    inSequences.buildEdgeGraph();
    std::vector<std::vector<Edge>> &adjEdgeList = inSequences.getAdjEdgeList();
     // get the headers to uIds table to look for the header
    
    phmap::flat_hash_map<std::string, unsigned int> &headersToIds = *inSequences.getHash1();
	std::vector<Path> alignmentPaths = inAlignments.getPaths(headersToIds);
    
    // map node names to internal ids
    NodeTable nodeTable(nodeFile, headersToIds);
    nodeTable.add(source, {headersToIds[source],1});
    nodeTable.add(destination, {headersToIds[destination],1});
    dist[pId] = 0;
    Path firstPath(nodeTable);
    firstPath.push_back(nodeTable[source].uId,'0');
    std::pair<const uint32_t,Path> *u = new std::pair<const uint32_t,Path>(pId++, firstPath);
    Q.insert(u, 0); // append source node
    lg.verbose("Starting search");
    while (Q.size() > 0 && steps < maxSteps) { // the main loop
        std::pair<const uint32_t,Path> u = *Q.extractMin(); // remove and return best segment
        bool pathFound = false;
        pathFound = evaluatePath(u.second, inSequences, alignmentPaths, dist, nodeTable, pId, nodeTable[destination].uId);
        if (pathFound)
            continue;
        
        uint32_t alt = dist[u.first] + 1;
        for(auto v : adjEdgeList.at(u.second.path.back().id)) {
            
//            if (v.id == nodeTable[destination].uId)
//                alt += 2000;
            
            if (u.second.path.back().orientation != '0' && u.second.path.back().orientation != v.orientation0)
                continue;
            
            if (u.second.path.back().orientation == '0') // set orientation of the start node for this path
                u.second.path.back().orientation = v.orientation0;
            
            InSegment &nextSegment = inSequences.findSegmentBySUId(v.id);
            lg.verbose("Inspecting segment: " + nextSegment.getSeqHeader() + v.orientation1);
            auto got = u.second.nodeTable.records.find(nextSegment.getSeqHeader());
            lg.verbose("We can visit this node n time: " + std::to_string(got->second.count));
            if (got != u.second.nodeTable.records.end() && got->second.count > 0) {
                Path newPath(u.second);
                newPath.push_back(v.id,v.orientation1);
                auto got2 = newPath.nodeTable.records.find(nextSegment.getSeqHeader());
                --got2->second.count;
                std::pair<const uint32_t,Path> *u = new std::pair<const uint32_t,Path>(pId++, newPath);
                Q.insert(u, alt);
                dist[pId] = alt;
            }
        }
        ++steps;
    }
    if (steps >= maxSteps)
        std::cout<<"Reached maximum number of steps ("<<+steps<<")"<<std::endl;
    lg.verbose("Search completed");
}
