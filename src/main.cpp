#include "main.h"

std::string version = "0.1";

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

int main(int argc, char **argv) {
    
    short int c; // optarg
    bool arguments = true;
    bool isPipe = false; // to check if input is from pipe
    unsigned long long int gSize = 0; // expected genome size, with 0 NG/LG* statistics are not computed
    std::string action, aligner = "GraphAligner", preset_type, preset = " -x vg", cmd;
    const char strHelp[] = "gfalign [options] [tool] [arguments]\n";
    
    if (argc == 1) { // gfastats with no arguments
        printf("%s", strHelp);
        printf("-h for additional help.\n");
        exit(0);
    }
    static struct option long_options[] = { // struct mapping long options
        {"input-sequence", required_argument, 0, 'f'},
        {"input-reads", required_argument, 0, 'r'},
        {"input-alignment", required_argument, 0, 'g'},
        {"node-list", required_argument, 0, 'n'},
        {"max-steps", required_argument, 0, 'm'},
        {"source", required_argument, 0, 's'},
        {"destination", required_argument, 0, 'd'},
        {"cmd", no_argument, &userInput.cmd_flag, 1},
        {"preset", required_argument, 0, 'p'},
        {"out-format", required_argument, 0, 'o'},
        {"sort-alignment", no_argument, &userInput.sortAlignment_flag, 1},
        
        {"output-terminal-alignments", no_argument, &userInput.terminalAlignments_flag, 1},

        {"threads", required_argument, 0, 'j'},
        
        {"verbose", no_argument, &verbose_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        
        {0, 0, 0, 0}
    };
    const static std::unordered_map<std::string,int> tools{
        {"align",1},
        {"eval",2},
        {"subgraph",3},
        {"dijkstra",4}
    };
    while (arguments) { // loop through argv
        
        int option_index = 0;
        c = getopt_long(argc, argv, "-:d:s:v:f:p:g:m:n:j:o:r:h",
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
                    case 2:
                    case 3:
                    case 4:
                        cmd = getArgs(optarg, argc, argv);
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
            case 'p': // presets
            {
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
            case 'g': // input alignment
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                    userInput.pipeType = 'g'; // pipe input is a sequence
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.inAlign = optarg;
                    userInput.stats_flag = true;
                }
                break;
            case 'm':
                userInput.dijkstraSteps = atoi(optarg);
                break;
            case 'n': // node list
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                    userInput.pipeType = 'g'; // pipe input is a sequence
                }else{ // input is a regular file
                    ifFileExists(optarg);
                    userInput.nodeList = optarg;
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
            case 's':
                userInput.source = optarg;
                break;
            case 'o': // handle output (file or stdout)
                userInput.outFile = optarg;
                break;
            case 'v': // software version
                printf("gfalign v%s\n", version.c_str());
                printf("Giulio Formenti giulio.formenti@gmail.com\n");
                exit(0);
                
            case 'h': // help
                printf("%s", strHelp);
                printf("\nTools:\n");
                printf("align\n");
                printf("eval\n");
                printf("\nOptions:\n");
                printf("-p --preset alignment presets (currently supports: hifi|CLR).\n");
                printf("-v --version software version.\n");
                printf("--cmd print $0 to stdout.\n");
                exit(0);
        }
        if (userInput.sortAlignment_flag || userInput.terminalAlignments_flag) // handle various cases in which the output should not include summary stats
            userInput.stats_flag = false;
    }
    if (userInput.cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++)
            printf("%s ", argv[arg_counter]);
        printf("\n");
    }
    if (tools.at(action) == 2){
        if (userInput.inAlign == "") {
            printf("%s", strHelp);
            printf("\nOptions:\n");
            printf("-f --input-sequence sequence input file (gfa1/2).\n");
            printf("-g --input-alignment alignment input file (currently supports: GAF).\n");
            printf("-o --out-format ouput to file or stdout (currently supports: GFA, GAF).\n");
            printf("--sort-alignment output sorted alignment.\n");
            printf("--output-terminal-alignments output terminal alignments.\n");
            exit(0);
        }
        Input in;
        in.load(userInput); // load user input
        lg.verbose("User input loaded");
        threadPool.init(maxThreads); // initialize threadpool
        lg.verbose("GFA: " + userInput.inSequence);
        InSequences inSequences; // initialize sequence collection object
        InAlignments inAlignments; // initialize alignment collection object
        lg.verbose("Alignment object generated");
        
        if(userInput.inAlign != ""){
            
            lg.verbose("Sequence object generated");
            in.read(inSequences); // read input content to inSequences container
            if (userInput.stats_flag) {
                Report report;
                report.reportStats(inSequences, gSize, 0);
            }
        }
        if(userInput.inAlign != ""){
            
            lg.verbose("Alignment: " + userInput.inAlign);
            in.read(inAlignments); // read input content to inAlignments container
            jobWait(threadPool);
            inAlignments.sortAlignmentsByNameAscending();
            inAlignments.markDuplicates();
            
            if(userInput.stats_flag)
                inAlignments.printStats();
            else if (userInput.sortAlignment_flag)
                inAlignments.outAlignments();
        }
        threadPool.join();
        if(userInput.inAlign != "" && userInput.outFile != ""){
            evalGFA(inSequences, inAlignments);
            Report report;
            report.writeToStream(inSequences, userInput.outFile, userInput);
        }
    }else if (tools.at(action) == 3){
        if (userInput.nodeList == "") {
            printf("%s", strHelp);
            printf("\nOptions:\n");
            printf("-f --input-sequence sequence input file (gfa1/2).\n");
            printf("-n --node-list list of nodes to retain in the subgraph.\n");
            printf("-o --out-format ouput to file or stdout (currently supports: GFA).\n");
            exit(0);
        }
        Input in;
        in.load(userInput); // load user input
        lg.verbose("User input loaded");
        threadPool.init(maxThreads); // initialize threadpool
        lg.verbose("GFA: " + userInput.inSequence);
        InSequences inSequences; // initialize sequence collection object
        lg.verbose("Sequence object generated");
        
        if(userInput.inSequence != ""){
            
            in.read(inSequences); // read input content to inSequences container
            std::vector<std::string> nodeList;
            std::string line; // Replace with your file's name
            std::ifstream file(userInput.nodeList);
            while (std::getline(file, line))
                nodeList.push_back(line);
            file.close();
            
            lg.verbose("Node list read");
            InSequences *subgraph = inSequences.subgraph(nodeList);
            Report report;
            if (userInput.outFile != "") // output sequences to file or stdout
                report.writeToStream(*subgraph, userInput.outFile, userInput);
            
            if (userInput.stats_flag) {
                Report report;
                report.reportStats(*subgraph, gSize, 0);
            }
            delete subgraph;
        }
        threadPool.join();
    }else if (tools.at(action) == 4){
        if (userInput.nodeList == "") {
            printf("%s", strHelp);
            printf("\nOptions:\n");
            printf("-f --input-sequence <filename> sequence input file (gfa1/2).\n");
            printf("-n --node-list <filename> list of nodes available to the search.\n");
            printf("-s --source <string> source node.\n");
            printf("-d --destination <string> destination node.\n");
            printf("-m --max-steps <int> limit graph exploration.\n");
            exit(0);
        }
        Input in;
        in.load(userInput); // load user input
        lg.verbose("User input loaded");
        threadPool.init(maxThreads); // initialize threadpool
        lg.verbose("GFA: " + userInput.inSequence);
        InSequences inSequences; // initialize sequence collection object
        lg.verbose("Sequence object generated");
        
        if(userInput.inSequence != ""){
            
            in.read(inSequences); // read input content to inSequences container
            std::vector<std::string> nodeList;
            std::string line; // Replace with your file's name
            std::ifstream file(userInput.nodeList);
            while (std::getline(file, line))
                nodeList.push_back(line);
            file.close();
            
            lg.verbose("Node list read");
            dijkstra(inSequences, nodeList, userInput.source, userInput.destination, userInput.dijkstraSteps);
        }
        threadPool.join();
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
