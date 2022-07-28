#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "functions.h" // global functions

#include "log.h"
#include "global.h"
#include "uid-generator.h"

#include "gfa-lines.h"
#include "gfa.h"

#include "alignments.h"

#include "eval.h"


void evalGFA(InSequences& InSequences, InAlignments& InAlignments) {
    
    InAlignments.buildEdgeGraph(InSequences.getHash1(), InSequences.getHash2());

    std::vector<std::vector<Edge>> adjEdgeList = InAlignments.getEdgeGraph();
    
    std::vector<InEdge>* edges = InSequences.getEdges();
    
    phmap::flat_hash_map<unsigned int, std::string>* idsToHeaders = InSequences.getHash2();
    
    for(InEdge& edge : *edges) {
        
        lg.verbose("Assessing edge: " + (*idsToHeaders)[edge.getsId1()] + "(" + std::to_string(edge.getsId1()) + ") " + edge.getsId1Or() + " " + (*idsToHeaders)[edge.getsId2()] + "(" + std::to_string(edge.getsId2()) + ") " + edge.getsId2Or());
        
        Edge fwEdge {edge.getsId1Or(), edge.getsId2(), edge.getsId2Or()};
        
        auto it = find(adjEdgeList.at(edge.getsId1()).begin(), adjEdgeList.at(edge.getsId1()).end(), fwEdge);
        
        if (it != adjEdgeList.at(edge.getsId1()).end()) {
            
            lg.verbose("Edge implied by read alignment (weight: " + std::to_string(it->weight) + ")");
            
            edge.appendTag({'i', "RC", std::to_string(it->weight)});
            
        }
        
    }

}
