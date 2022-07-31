#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <unordered_map>
#include <getopt.h>
#include <vector>

#include <iostream>
#include <fstream>

#include <parallel_hashmap/phmap.h>


#include "bed.h"

#include "struct.h"
#include "functions.h"
#include "global.h"

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "struct.h"
#include "stream-obj.h"

#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"

#include "alignments.h"
#include "eval.h"
#include "input.h"
#include "output.h"

#include <main.h>

std::string version = "0.1";

//global
std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

int hc_flag;
int hc_cutoff;
short int tabular_flag;
int cmd_flag;
int verbose_flag;
int outBubbles_flag;
int stats_flag;
int discoverPaths_flag;
int outFile_flag;
int maxThreads = 0;

std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
Log lg;


int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    bool isPipe = false; // to check if input is from pipe
    
    unsigned long long int gSize = 0; // expected genome size, with 0 NG/LG* statistics are not computed
    
    UserInput userInput; // initialize input object
    
    std::string action, aligner, cmd;
    
    userInput.outSequence = "gfa"; // default output type
    
    if (argc == 1) { // gfastats with no arguments
            
        printf("gfalign [command]\n-h for additional help.\n");
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"input-sequence", required_argument, 0, 'f'},
        {"input-reads", required_argument, 0, 'r'},
        {"input-alignment", required_argument, 0, 'g'},
        {"cmd", no_argument, &cmd_flag, 1},
        {"aligner", required_argument, 0, 'a'},
        {"out-format", required_argument, 0, 'o'},
        {"no-sequence", no_argument, &userInput.noSequence, 1},

        {"threads", required_argument, 0, 'j'},
        
        {"verbose", no_argument, &verbose_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        
        {0, 0, 0, 0}
    };
    
    const static std::unordered_map<std::string,int> tools{
        {"align",1},
        {"eval",2}
    };
    
    while (arguments) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:v:f:a:g:j:o:r:h",
                        long_options, &option_index);
        
        if (optind < argc && !isPipe) { // if pipe wasn't assigned already
            
            isPipe = isDash(argv[optind]) ? true : false; // check if the argument to the option is a '-' and set it as pipe input
            
        }
        
        if (optarg != nullptr && !isPipe) { // case where pipe input is given as positional argument (input sequence file)
        
            isPipe = isDash(optarg) ? true : false;
            
        }
        
        if (c == -1) { // exit the loop if run out of options
            break;
            
        }

        switch (c) {
            case ':': // handle options without arguments
                switch (optopt) { // the command line option last matched
                    case 'a':
                        aligner = "GraphAligner";
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
                        
                        cmd = "GraphAligner" + getArgs(optarg, argc, argv);
                        
                        std::cout<<"Invoking: "<<cmd<<std::endl;
                        std::system(cmd.c_str());
                        
                        arguments = false;
                        
                        break;
                    case 2:
                        
                        cmd = getArgs(optarg, argc, argv);
                        
                        break;
                        
                }
                
            case 0: // case for long options without short options
                
//                if (strcmp(long_options[option_index].name,"line-length") == 0)
//                  splitLength = atoi(optarg);
                
                break;
                
            case 'f': // input sequence
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'f'; // pipe input is a sequence
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.iSeqFileArg = optarg;
                    stats_flag = true;
                    
                }
                    
                break;

            case 'g': // input alignment
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'g'; // pipe input is a sequence
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.iAlignFileArg = optarg;
                    stats_flag = true;
                    
                }
                    
                break;
                
            case 'j': // max threads
                maxThreads = atoi(optarg);
                stats_flag = 1;
                break;
                
            case 'r': // input reads
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'r'; // pipe input is a sequence
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.iReadFileArg = optarg;
                    stats_flag = true;
                    
                }
                    
                break;
            case 'o': // handle output (file or stdout)
                userInput.outSequence = optarg;
                outFile_flag = 1;
                break;
                
            case 'v': // software version
                printf("gfalign v%s\n", version.c_str());
                printf("Giulio Formenti giulio.formenti@gmail.com\n");
                exit(0);
                
            case 'h': // help
                printf("gfalign [command]\n");
                printf("\nOptions:\n");
                printf("-v --version software version.\n");
                printf("--cmd print $0 to stdout.\n");
                exit(0);
        }
        
        if    (argc == 2 || // handle various cases in which the output should include summary stats
              (argc == 3 && pos_op == 2) ||
              (argc == 4 && pos_op == 3)) {
            
        }
        
    }
    
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
    
    if (tools.at(action) == 2){
        
        StreamObj streamObj;
        std::shared_ptr<std::istream> stream;
        std::string newLine;
        
        Input in;
        
        in.load(userInput); // load user input
        
        lg.verbose("User input loaded");
        
        threadPool.init(maxThreads); // initialize threadpool
        
        lg.verbose("GFA: " + userInput.iSeqFileArg);
    
        InSequences inSequences; // initialize sequence collection object
        
        InAlignments inAlignments; // initialize alignment collection object
        
        lg.verbose("Alignment object generated");
        
        if(userInput.iSeqFileArg != ""){
            
            lg.verbose("Sequence object generated");
            
            in.read(inSequences); // read input content to inSequences container
            
            jobWait(threadPool);
            
            Report report;
            
            report.reportStats(inSequences, gSize);
        
        }
        
        if(userInput.iAlignFileArg != ""){
            
            lg.verbose("Alignment: " + userInput.iAlignFileArg);
            
            in.read(inAlignments); // read input content to inAlignments container
            
            jobWait(threadPool);
            
            inAlignments.printStats();
            
        }
        
        threadPool.join();
        
        if(userInput.iSeqFileArg != "" && userInput.iAlignFileArg != ""){
            
            evalGFA(inSequences, inAlignments);
            
            Report report;
            
            if (outFile_flag) { // output sequences to file or stdout
            
                report.outFile(inSequences, userInput);
                
            }
            
        }
        
    }
    
    exit(EXIT_SUCCESS);
    
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
