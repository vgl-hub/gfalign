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

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "functions.h" // global functions

#include "log.h"
#include "global.h"
#include "uid-generator.h"

#include "gfa-lines.h"
#include "reads.h"

#include "threadpool.h"
#include "gfa.h"
#include "sak.h" // swiss army knife

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "stream-obj.h"

#include "alignments.h"

#include "input.h"

void Input::load(UserInput userInput) {
    
    this->userInput = userInput;
    
}

void Input::read(InAlignments& inAlignments) {
    
    if (userInput.iAlignFileArg.empty()) {return;}
        
    threadPool.init(maxThreads); // initialize threadpool

    inAlignments.load(userInput);
    
    while (true) {
        
        if (threadPool.empty()) {threadPool.join(); break;}
        lg.verbose("Remaining jobs: " + std::to_string(threadPool.queueSize()), true);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        
    }

}

void Input::read(InSequences& inSequences) {
    
    if (userInput.iAlignFileArg.empty()) {return;}
        
    threadPool.init(maxThreads); // initialize threadpool
    
    readGFA(inSequences, userInput, stream);
    
    while (true) {
        
        if (threadPool.empty()) {threadPool.join(); break;}
        lg.verbose("Remaining jobs: " + std::to_string(threadPool.queueSize()), true);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        
    }

}
