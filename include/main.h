#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <unistd.h>
#include <getopt.h>

#include <vector>  //required for zstream
#include <stack>
#include <queue>
#include <string.h>
#include <algorithm> //required for zstream
#include <cstring> //required for zstream
#include <tuple> // for graph manipulation
#include <cctype> // toupper()
#include <iomanip>
#include <map>

#include <chrono>
#include <memory>

#include <thread>
#include <mutex>
#include <condition_variable>

#include "log.h"

#include "uid-generator.h"

#include "bed.h"

#include "global.h" // global variables
#include "struct.h"
#include "functions.h" // global functions

#include "threadpool.h"

#include <parallel-hashmap/phmap.h>

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/ozstream.hpp>
#include <zstream/ozstream_impl.hpp>

#include "gfa-lines.h"

#include "gfa.h" // gfa classes
#include "sak.h" // swiss army knife classes

#include "stream-obj.h"

#include "output.h" // output classes
#include "alignments.h"
#include "input.h"
#include "eval.h"

std::string getArgs(char* optarg, unsigned int argc, char **argv);

#endif /* MAIN_H */
