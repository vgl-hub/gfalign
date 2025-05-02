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

struct PathAlignmentStats {
	uint32_t badAlignments = 0, goodAlignments = 0, unaligned = 0;
};

PathAlignmentStats evaluatePath(const Path &path, InSequences &inSequences, std::vector<Path> alignmentPaths, bool filterAlignments = false, bool printAlignments = false) {
	
	Step lastStep = path.path.back();
	InSegment &segment = inSequences.findSegmentBySUId(lastStep.id);
	lg.verbose("We are at segment: " + segment.getSeqHeader() + lastStep.orientation);
	if (printAlignments)
		std::cout<<path.returnPath(inSequences)<<std::endl;

	PathAlignmentStats pathAlignmentStats;
	phmap::flat_hash_set<uint32_t> uIds;
	for (Step step : path.path)
		uIds.insert(step.id);
	int dp[MAX_N][MAX_N] = {{0}};
	for (Path &alignmentPath : alignmentPaths) {
		if (filterAlignments) {
			bool next = false;
			for (Step step : alignmentPath.path) { // remove spurious reads
				if (uIds.find(step.id) == uIds.end()) {
					next = true;
					++pathAlignmentStats.unaligned;
				}
			}
			if (next)
				continue;
		}
		PairwisePathAlignment alignmentFw = alignPaths(0, -1, -1, path, alignmentPath, dp);
		PairwisePathAlignment alignmentRc = alignPaths(0, -1, -1, path, alignmentPath.reverseComplement(), dp);
		int32_t bestAlignmentScore = (alignmentFw.alignmentScore > alignmentRc.alignmentScore) ? alignmentFw.alignmentScore : alignmentRc.alignmentScore;
		if (bestAlignmentScore < 0)
			++pathAlignmentStats.badAlignments;
		else
			++pathAlignmentStats.goodAlignments;
		
		if (printAlignments) {
			std::string aln = (alignmentFw.alignmentScore > alignmentRc.alignmentScore) ? alignmentFw.getAlignment(*inSequences.getHash2(), true) : alignmentRc.getAlignment(*inSequences.getHash2(), true);
			std::cout<<aln<<'\t'<<alignmentPath.qName<<'\t'<<bestAlignmentScore<<std::endl;
//			std::cout<<alignmentFw.getAlignment(*inSequences.getHash2())<<" "<<alignmentFw.alignmentScore<<std::endl;
//			std::cout<<alignmentRc.getAlignment(*inSequences.getHash2())<<" "<<alignmentRc.alignmentScore<<std::endl;
		}
	}
	return pathAlignmentStats;
}

void dijkstra(InSequences &inSequences, InAlignments& inAlignments, std::string nodeFile, std::string source, std::string destination, uint32_t maxSteps, uint32_t minNodes, bool returnAllPaths = false) {
	
	Path bestPath;
	int32_t bestPath_alt = std::numeric_limits<int>::max();
	uint64_t pathCounter = 0;
	uint32_t bestPath_uniques = 0, steps = 0, pId = 0;
	std::vector<uint64_t> destinations;
	FibonacciHeap<std::pair<const uint32_t,Path>*> Q; // node priority queue Q
	inSequences.buildEdgeGraph();
	std::vector<std::vector<Edge>> &adjEdgeList = inSequences.getAdjEdgeList();
	 // get the headers to uIds table to look for the header
	
	phmap::flat_hash_map<std::string, unsigned int> &headersToIds = *inSequences.getHash1();
	std::vector<Path> alignmentPaths = inAlignments.getPaths(headersToIds);
	
	// map node names to internal ids
	NodeTable nodeTable(nodeFile, headersToIds);
	nodeTable.add(source, {headersToIds[source],1});
	nodeTable.add(destination, {headersToIds[destination],1});
	Path firstPath(nodeTable);
	firstPath.push_back(nodeTable[source].uId,'0');
	std::pair<const uint32_t,Path> *u = new std::pair<const uint32_t,Path>(pId++, firstPath);
	Q.insert(u, 0); // append source node
	lg.verbose("Starting search");
	while (Q.size() > 0 && steps < maxSteps) { // the main loop
		std::pair<const uint32_t,Path> *u = Q.extractMin(); // remove and return best path to extend
		for(auto v : adjEdgeList.at(u->second.path.back().id)) {
			if (u->second.path.back().orientation != '0' && u->second.path.back().orientation != v.orientation0)
				continue;
			
			InSegment &nextSegment = inSequences.findSegmentBySUId(v.id);
			lg.verbose("Inspecting segment: " + nextSegment.getSeqHeader() + v.orientation1);
			auto got = u->second.nodeTable.records.find(nextSegment.getSeqHeader());
			if (got != u->second.nodeTable.records.end() && got->second.count > 0) {
				
				lg.verbose("We can visit this node n times: " + std::to_string(got->second.count));
				Path newPath(u->second);
			
				if (newPath.path.back().orientation == '0') // set orientation of the start node for this path
					newPath.path.back().orientation = v.orientation0;
				
				newPath.push_back(v.id,v.orientation1);
				
				std::vector<std::string> uniques;
				for (auto n : newPath.path) {
					InSegment &segment = inSequences.findSegmentBySUId(n.id);
					uniques.push_back(segment.getSeqHeader());
				}
				std::sort(uniques.begin(), uniques.end());
				auto last = std::unique(uniques.begin(), uniques.end());
				uniques.erase(last, uniques.end());
				
				PathAlignmentStats pathAlignmentStats = evaluatePath(newPath, inSequences, alignmentPaths, true);
				int32_t alt = (int32_t)pathAlignmentStats.badAlignments - (int32_t)pathAlignmentStats.goodAlignments - (int32_t)uniques.size();
				
				if(v.id != nodeTable[destination].uId) {
					auto got2 = newPath.nodeTable.records.find(nextSegment.getSeqHeader());
					--got2->second.count;
					std::pair<const uint32_t,Path> *u = new std::pair<const uint32_t,Path>(pId++, newPath);
					Q.insert(u, alt);
				}else{
					lg.verbose("Destination found.");
					++pathCounter;
					phmap::flat_hash_map<uint32_t,uint32_t> pathNodes = newPath.pathToMap();
					bool hamiltonian = nodeTable.checkHamiltonian(pathNodes, newPath.size()); // check if hamiltonian
					bool printPath = false;
					if (uniques.size() >= minNodes && (bestPath_uniques < uniques.size() || (bestPath_uniques == uniques.size() && bestPath_alt > alt))) { // better path found
						bestPath = newPath;
						bestPath_alt = alt;
						bestPath_uniques = uniques.size();
						printPath = true;
					}
					if (returnAllPaths || printPath)
						std::cout<<+pathCounter<<'\t'<<+pathAlignmentStats.badAlignments<<'\t'<<+pathAlignmentStats.goodAlignments<<'\t'<<+alt<<'\t'<<newPath.size()<<'\t'<<uniques.size()<<'\t'<<(hamiltonian ? 'T' : 'F')<<'\t'<<newPath.returnPath(inSequences)<<std::endl; // print new path only if better
				}
			}
		}
		++steps;
		delete u;
	}
	if (steps >= maxSteps)
		std::cout<<"Reached maximum number of steps ("<<+steps<<")"<<std::endl;
	lg.verbose("Search completed");
}


void evalPath(InSequences &inSequences, InAlignments& inAlignments, std::string pathStr) {
	Path path;
	std::vector<std::string> components;
	char sIdOr;
	
	phmap::flat_hash_map<std::string, unsigned int> &headersToIds = *inSequences.getHash1();
	inSequences.uId.next();
	std::vector<char> delimiters {';', ','};
	components = readDelimitedArr(pathStr, delimiters, "", true);
	
	for (auto it = std::begin(components); it != std::end(components); ++it) {
		
		if(it == std::begin(components) && *it == "") { // handle starting/ending gap
			fprintf(stderr, "Error: cannot handle starting gap. Terminating.\n");
			exit(1);
		}
		
		std::string component = *it;
		if (std::next(it) != std::end(components))
			component.pop_back(); // remove separator
		sIdOr = component.back(); // get sequence orientation
		component.pop_back();
				
		auto got = headersToIds.find(component); // get the headers to uIds table (remove sequence orientation in the gap first)
	
		if (got == headersToIds.end()) { // this is the first time we see this segment
			fprintf(stderr, "Error: cannot find node (%s). Terminating.\n", component.c_str());
			exit(1);
		}else{
			path.push_back(got->second,sIdOr);
		}
	}
	std::vector<std::string> uniques;
	for (auto n : path.path) {
		InSegment &segment = inSequences.findSegmentBySUId(n.id);
		uniques.push_back(segment.getSeqHeader());
	}
	std::sort(uniques.begin(), uniques.end());
	auto last = std::unique(uniques.begin(), uniques.end());
	uniques.erase(last, uniques.end());
	
	std::vector<Path> alignmentPaths = inAlignments.getPaths(headersToIds);
	PathAlignmentStats pathAlignmentStats = evaluatePath(path, inSequences, alignmentPaths, false, true);
	int32_t alt = (int32_t)pathAlignmentStats.badAlignments - (int32_t)pathAlignmentStats.goodAlignments - (int32_t)uniques.size();
	
	std::cout<<+pathAlignmentStats.badAlignments<<"\t"<<+pathAlignmentStats.goodAlignments<<"\t"<<+alt<<"\t"<<path.size()<<"\t"<<uniques.size()<<std::endl;
}
