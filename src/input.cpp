#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>

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

#include "alignments.h"
#include "input.h"

void Input::load(UserInputGfalign userInput) {
    this->userInput = userInput;
}

void Input::read(InAlignments& inAlignments) {
    
    if (userInput.inAlign.empty()) {return;}
    
    std::shared_ptr<std::istream> stream;
    stream = streamObj.openStream(userInput, 'g');
    inAlignments.load(stream, userInput.terminalAlignments_flag);
}

void Input::read(InSequences& inSequences) {
    
    if (userInput.inAlign.empty()) {return;}
    
    stream = streamObj.openStream(userInput, 'f');
    readGFA(inSequences, userInput, stream);
    jobWait(threadPool);
    inSequences.updateStats();
}
