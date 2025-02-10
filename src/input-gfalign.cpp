#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>

#include <iostream>
#include <fstream>
#include <sstream>

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

#include "nodetable.h"
#include "alignments.h"
#include "eval.h"
#include "output.h" // output classes
#include "input-gfalign.h"

void Input::loadInput(UserInputGfalign userInput) {
    this->userInput = userInput;
	
	if(userInput.inSequence != ""){
		lg.verbose("GFA: " + userInput.inSequence);
		stream = streamObj.openStream(this->userInput, 'f');
		readGFA(this->inSequences, userInput, stream);
		jobWait(threadPool);
		inSequences.updateStats();
		lg.verbose("Sequence object generated");
		if (userInput.stats_flag) {
			Report report;
			report.reportStats(inSequences, 0, 0);
		}
	}
	if(userInput.inAlign != ""){
		lg.verbose("Alignment: " + userInput.inAlign);
		read(this->inAlignments); // read input content to inAlignments container
		jobWait(threadPool);
	}
}

void Input::read(InAlignments& inAlignments) {
    
    if (userInput.inAlign.empty()) {return;}
    inAlignments.load(userInput.inAlign, userInput.terminalAlignments_flag);
}

void Input::read() {
    
    switch (userInput.mode) {
        case 0: { // graph alignment
            
        }
        case 1: { // graph evaluation
            lg.verbose("GFA: " + userInput.inSequence);
            InSequences inSequences; // initialize sequence collection object
            InAlignments inAlignments; // initialize alignment collection object
            lg.verbose("Alignment object generated");

            if(userInput.inAlign != ""){
                
                lg.verbose("Alignment: " + userInput.inAlign);
                read(inAlignments); // read input content to inAlignments container
                jobWait(threadPool);
                inAlignments.sortAlignmentsByNameAscending();
                inAlignments.markDuplicates();
                
                if(userInput.alignStats_flag)
                    inAlignments.printStats();
                else if (userInput.sortAlignment_flag)
                    inAlignments.outAlignments();
            }
            if(userInput.inAlign != "" && userInput.outFile != ""){
                evalGFA(this->inSequences, inAlignments);
                Report report;
                report.writeToStream(this->inSequences, userInput.outFile, userInput);
            }
            break;
        }
        case 2: { // subgraph

            std::vector<std::string> nodeList;
            std::string line; // Replace with your file's name
            std::ifstream file(userInput.nodeFile);
            while (std::getline(file, line))
                nodeList.push_back(line);
            file.close();
            
            lg.verbose("Node list read");
            InSequences *subgraph = inSequences.subgraph(nodeList);
            Report report;
            if (userInput.outFile != "") // output sequences to file or stdout
                report.writeToStream(*subgraph, userInput.outFile, userInput);
            delete subgraph;
        break;
        }
        case 3: { // path search
			dijkstra(this->inSequences, this->inAlignments, userInput.nodeFile, userInput.source, userInput.destination, userInput.dijkstraSteps);
            break;
        }
    }
}
