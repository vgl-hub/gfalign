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
#include "gfa-lines.h"
#include "gfa.h" // gfa classes
#include "sak.h" // swiss army knife classes
#include "stream-obj.h"

#include <parallel-hashmap/phmap.h>
#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/ozstream.hpp>
#include <zstream/ozstream_impl.hpp>

#include "alignments.h"
#include "input-gfalign.h"
#include "main.h"

std::string version = "0.0.1";

//global
std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

int verbose_flag;
Log lg;
std::vector<Log> logs;
int tabular_flag;
int maxThreads = 0;
std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
UserInputGfalign userInput; // initialize input object

void printHelp() {
    printf("gfalign [options] [tool] [arguments]\n-h for additional help.\n");
    printf("\nTools:\n");
    printf("align\n");
    printf("eval\n");
    printf("subgraph\n");
    printf("search\n");
    exit(0);
}

std::string getArgs(char* optarg, unsigned int argc, char **argv) {
    
    std::string cmd;
    bool record = false;

    for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
        
        if (optarg != argv[arg_counter] && !record) {
            continue;
        }else{
            record = true;
            if(optarg != argv[arg_counter]){
                cmd += ' ';
                cmd += argv[arg_counter];
            }
        }
    }
    return cmd;
}

int main(int argc, char **argv) {
    
    short int c; // optarg
    bool arguments = true;
    bool isPipe = false; // to check if input is from pipe

    if (argc == 1) // gfalign with no arguments
        printHelp();
    const static std::unordered_map<std::string,int> tools{
        {"align",0},
        {"eval",1},
        {"subgraph",2},
        {"search",3}
    };
    
    auto got = tools.find(argv[1]);
    if (got != tools.end()) {
        userInput.mode = got->second;
    }else{
        fprintf(stderr, "mode '%s' does not exist. Terminating\n", argv[1]);
        return EXIT_FAILURE;
    }
    
    switch (userInput.mode) {
        case 0: { // graph alignment
            std::string action, aligner = "GraphAligner", preset_type, preset = " -x vg", cmd;
            static struct option long_options[] = { // struct mapping long options
                {"input-sequence", required_argument, 0, 'f'},
                {"input-alignment", required_argument, 0, 'g'},
                {"preset", required_argument, 0, 'p'},
                {"input-reads", required_argument, 0, 'r'},
                {"out-format", required_argument, 0, 'o'},
                
                {"sort-alignment", no_argument, &userInput.sortAlignment_flag, 1},
                {"output-terminal-alignments", no_argument, &userInput.terminalAlignments_flag, 1},
                
                {"threads", required_argument, 0, 'j'},
                {"cmd", no_argument, &userInput.cmd_flag, 1},
                {"verbose", no_argument, &verbose_flag, 1},
                {"version", no_argument, 0, 'v'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };
            while (arguments) { // loop through argv
                
                int option_index = 1;
                c = getopt_long(argc, argv, "-:f:g:h:j:o:p:r:v",
                                long_options, &option_index);
                
                if (optind < argc && !isPipe) // if pipe wasn't assigned already
                    isPipe = isDash(argv[optind]) ? true : false; // check if the argument to the option is a '-' and set it as pipe input
                
                if (optarg != nullptr && !isPipe) // case where pipe input is given as positional argument (input sequence file)
                    isPipe = isDash(optarg) ? true : false;
                
                if (c == -1) // exit the loop if run out of options
                    break;
                
                switch (c) {
                    case ':': // handle options without arguments
                        switch (optopt) { // the command line option last matched
                            case 'p':
                                aligner = "GraphAligner";
                                preset = " --seeds-mxm-length 1000 --min-alignment-score 1000 --precise-clipping 0.75";
                                break;
                            default:
                                fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                                return EXIT_FAILURE;
                        }
                        break;
                    default: // handle positional arguments
                        
                        action = optarg;
                        switch (tools.count(optarg) ? tools.at(optarg) : 0) {
                            case 1:
                                cmd = aligner + getArgs(optarg, argc, argv) + preset;
                                std::cout<<"Invoking: "<<cmd<<std::endl;
                                std::system(cmd.c_str());
                                arguments = false;
                                break;
                            default:
                                std::cout<<"Could not find command: "<<optarg<<std::endl;
                                exit(1);
                                break;
                        }
                    case 0: // case for long options without short options
                        
                        //                if (strcmp(long_options[option_index].name,"line-length") == 0)
                        //                  splitLength = atoi(optarg);
                        break;
                    case 'p': { // presets
                    
                        const std::map<std::string, std::vector<std::string>> presets {
                            {"hifi", {"GraphAligner", " -x vg"}},
                            {"CLR", {"GraphAligner", " -x vg --seeds-mxm-length 1000 --min-alignment-score 1000 --precise-clipping 0.75"}}
                        };
                        if (presets.find(optarg) != presets.end()){
                            
                            aligner = presets.at(optarg)[0];
                            preset = presets.at(optarg)[1];
                        }else{
                            
                            std::cout<<"Could not find preset: "<<optarg<<std::endl;
                            exit(1);
                        }
                        break;
                    }
                    case 'f': // input sequence
                        
                        if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                            userInput.pipeType = 'f'; // pipe input is a sequence
                        }else{ // input is a regular file
                            
                            ifFileExists(optarg);
                            userInput.inSequence = optarg;
                            userInput.stats_flag = true;
                        }
                        break;
                    case 'g': // input alignment
                        if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                            userInput.pipeType = 'g'; // pipe input is a sequence
                        }else{ // input is a regular file
                            
                            ifFileExists(optarg);
                            userInput.inAlign = optarg;
                            userInput.stats_flag = true;
                        }
                        break;
                    case 'j': // max threads
                        maxThreads = atoi(optarg);
                        userInput.stats_flag = 1;
                        break;
                        
                    case 'r': // input reads
                        if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                            
                            userInput.pipeType = 'r'; // pipe input is a sequence
                        }else{ // input is a regular file
                            optind--;
                            for( ;optind < argc && *argv[optind] != '-'; optind++){
                                ifFileExists(argv[optind]);
                                userInput.inFiles.push_back(argv[optind]);
                            }
                        }
                        break;
                    case 'o': // handle output (file or stdout)
                        userInput.outFile = optarg;
                        break;
                    case 'v': // software version
                        printf("gfalign v%s\n", version.c_str());
                        printf("Giulio Formenti giulio.formenti@gmail.com\n");
                        exit(0);
                        
                    case 'h': // help
                        printf("gfalign align [options]\n");
                        printf("\nOptions:\n");
                        printf("-f --input-sequence sequence input file (gfa1/2).\n");
                        printf("-g --input-alignment alignment input file (currently supports: GAF).\n");
                        printf("-o --out-format ouput to file or stdout (currently supports: GAF).\n");
                        printf("-p --preset alignment presets (currently supports: hifi|CLR).\n");
                        printf("-v --version software version.\n");
                        printf("--cmd print $0 to stdout.\n");
                        exit(0);
                }
                if (userInput.sortAlignment_flag || userInput.terminalAlignments_flag) // handle various cases in which the output should not include summary stats
                    userInput.stats_flag = false;
            }
            break;
        }
        case 1: { // graph evaluation
            static struct option long_options[] = { // struct mapping long options
                {"input-sequence", required_argument, 0, 'f'},
                {"input-alignment", required_argument, 0, 'g'},
                {"out-format", required_argument, 0, 'o'},
                
                {"sort-alignment", no_argument, &userInput.sortAlignment_flag, 1},
                {"output-terminal-alignments", no_argument, &userInput.terminalAlignments_flag, 1},
                
                {"threads", required_argument, 0, 'j'},
                {"cmd", no_argument, &userInput.cmd_flag, 1},
                {"verbose", no_argument, &verbose_flag, 1},
                {"version", no_argument, 0, 'v'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };
            while (arguments) { // loop through argv
                
                int option_index = 1;
                c = getopt_long(argc, argv, "-:f:g:h:j:o:v",
                                long_options, &option_index);
                
                if (optind < argc && !isPipe) // if pipe wasn't assigned already
                    isPipe = isDash(argv[optind]) ? true : false; // check if the argument to the option is a '-' and set it as pipe input
                
                if (optarg != nullptr && !isPipe) // case where pipe input is given as positional argument (input sequence file)
                    isPipe = isDash(optarg) ? true : false;
                
                if (c == -1) // exit the loop if run out of options
                    break;
                
                switch (c) {
                    case ':': // handle options without arguments
                        break;
                    default: // handle positional arguments
                    case 0: // case for long options without short options
                        
                        //                if (strcmp(long_options[option_index].name,"line-length") == 0)
                        //                  splitLength = atoi(optarg);
                        break;
                    case 'f': // input sequence
                        
                        if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                            userInput.pipeType = 'f'; // pipe input is a sequence
                        }else{ // input is a regular file
                            ifFileExists(optarg);
                            userInput.inSequence = optarg;
                            userInput.stats_flag = true;
                        }
                        break;
                    case 'g': // input alignment
                        if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                            userInput.pipeType = 'g'; // pipe input is a sequence
                        }else{ // input is a regular file
                            
                            ifFileExists(optarg);
                            userInput.inAlign = optarg;
                            userInput.stats_flag = true;
                        }
                        break;
                        
                    case 'j': // max threads
                        maxThreads = atoi(optarg);
                        userInput.stats_flag = 1;
                        break;
                    case 'o': // handle output (file or stdout)
                        userInput.outFile = optarg;
                        break;
                    case 'v': // software version
                        printf("gfalign v%s\n", version.c_str());
                        printf("Giulio Formenti giulio.formenti@gmail.com\n");
                        exit(0);
                    case 'h': // help
                        printf("gfalign eval [options]\n");
                        printf("\nOptions:\n");
                        printf("-f --input-sequence sequence input file (gfa1/2).\n");
                        printf("-g --input-alignment alignment input file (currently supports: GAF).\n");
                        printf("-o --out-format ouput to file or stdout (currently supports: GFA, GAF).\n");
                        printf("--sort-alignment output sorted alignment.\n");
                        printf("--output-terminal-alignments output terminal alignments.\n");
                        exit(0);
                }
            }
            break;
        }
        case 2: { // subgraph
            static struct option long_options[] = { // struct mapping long options
                {"input-sequence", required_argument, 0, 'f'},
                {"node-file", required_argument, 0, 'n'},
                {"out-format", required_argument, 0, 'o'},
                
                {"threads", required_argument, 0, 'j'},
                {"cmd", no_argument, &userInput.cmd_flag, 1},
                {"verbose", no_argument, &verbose_flag, 1},
                {"version", no_argument, 0, 'v'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };
            while (arguments) { // loop through argv
                
                int option_index = 1;
                c = getopt_long(argc, argv, "-:f:h:j:n:o:v",
                                long_options, &option_index);
                
                if (optind < argc && !isPipe) // if pipe wasn't assigned already
                    isPipe = isDash(argv[optind]) ? true : false; // check if the argument to the option is a '-' and set it as pipe input
                
                if (optarg != nullptr && !isPipe) // case where pipe input is given as positional argument (input sequence file)
                    isPipe = isDash(optarg) ? true : false;
                
                if (c == -1) // exit the loop if run out of options
                    break;
                
                switch (c) {
                    case ':': // handle options without arguments
                        break;
                    default: // handle positional arguments
                    case 0: // case for long options without short options
                        
                        //                if (strcmp(long_options[option_index].name,"line-length") == 0)
                        //                  splitLength = atoi(optarg);
                        break;
                    case 'f': // input sequence
                        if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                            userInput.pipeType = 'f'; // pipe input is a sequence
                        }else{ // input is a regular file
                            
                            ifFileExists(optarg);
                            userInput.inSequence = optarg;
                            userInput.stats_flag = true;
                        }
                        break;
                    case 'j': // max threads
                        maxThreads = atoi(optarg);
                        userInput.stats_flag = 1;
                        break;
                    case 'o': // handle output (file or stdout)
                        userInput.outFile = optarg;
                        break;
                    case 'v': // software version
                        printf("gfalign v%s\n", version.c_str());
                        printf("Giulio Formenti giulio.formenti@gmail.com\n");
                        exit(0);
                    case 'h': // help
                        printf("gfalign subgraph [options]");
                        printf("\nOptions:\n");
                        printf("-f --input-sequence sequence input file (gfa1/2).\n");
                        printf("-n --node-file list of nodes to retain in the subgraph.\n");
                        printf("-o --out-format ouput to file or stdout (currently supports: GFA).\n");
                        exit(0);
                }
                break;
            }
        }
        case 3: { // path search
            static struct option long_options[] = { // struct mapping long options
                {"destination", required_argument, 0, 'd'},
                {"input-sequence", required_argument, 0, 'f'},
                {"max-steps", required_argument, 0, 'm'},
                {"node-file", required_argument, 0, 'n'},
                {"out-format", required_argument, 0, 'o'},
                {"source", required_argument, 0, 's'},
                
                {"threads", required_argument, 0, 'j'},
                {"cmd", no_argument, &userInput.cmd_flag, 1},
                {"verbose", no_argument, &verbose_flag, 1},
                {"version", no_argument, 0, 'v'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };
            while (arguments) { // loop through argv
                
                int option_index = 1;
                c = getopt_long(argc, argv, "-:f:h:j:n:o:v",
                                long_options, &option_index);
                
                if (optind < argc && !isPipe) // if pipe wasn't assigned already
                    isPipe = isDash(argv[optind]) ? true : false; // check if the argument to the option is a '-' and set it as pipe input
                
                if (optarg != nullptr && !isPipe) // case where pipe input is given as positional argument (input sequence file)
                    isPipe = isDash(optarg) ? true : false;
                
                if (c == -1) // exit the loop if run out of options
                    break;
                
                switch (c) {
                    case ':': // handle options without arguments
                        break;
                    default: // handle positional arguments
                    case 0: // case for long options without short options
                        
                        //                if (strcmp(long_options[option_index].name,"line-length") == 0)
                        //                  splitLength = atoi(optarg);
                        break;
                    case 'd':
                        userInput.destination = optarg;
                        break;
                    case 'f': // input sequence
                        if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                            userInput.pipeType = 'f'; // pipe input is a sequence
                        }else{ // input is a regular file
                            
                            ifFileExists(optarg);
                            userInput.inSequence = optarg;
                            userInput.stats_flag = true;
                        }
                        break;
                    case 'm':
                        userInput.dijkstraSteps = atoi(optarg);
                        break;
                    case 'n': // node list
                        ifFileExists(optarg);
                        userInput.nodeFile = optarg;
                        userInput.stats_flag = true;
                        break;
                    case 's':
                        userInput.source = optarg;
                        break;
                    case 'j': // max threads
                        maxThreads = atoi(optarg);
                        userInput.stats_flag = 1;
                        break;
                    case 'o': // handle output (file or stdout)
                        userInput.outFile = optarg;
                        break;
                    case 'v': // software version
                        printf("gfalign v%s\n", version.c_str());
                        printf("Giulio Formenti giulio.formenti@gmail.com\n");
                        exit(0);
                    case 'h': // help
                        printf("gfalign search [options]");
                        printf("\nOptions:\n");
                        printf("-f --input-sequence <filename> sequence input file (gfa1/2).\n");
                        printf("-n --node-file <filename> list of nodes available to the search.\n");
                        printf("-s --source <string> source node.\n");
                        printf("-d --destination <string> destination node.\n");
                        printf("-m --max-steps <int> limit graph exploration.\n");
                        exit(0);
                }
            }
            break;
        }
    }
    if (userInput.cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
    }
    Input in;
    threadPool.init(maxThreads); // initialize threadpool
    maxMem = (userInput.maxMem == 0 ? get_mem_total(3) * 0.9 : userInput.maxMem); // set memory limit
    in.loadInput(userInput); // load user input
    lg.verbose("User input loaded");
    in.read(); // read input reads and validate
    threadPool.join(); // join threads
    exit(EXIT_SUCCESS);
}
