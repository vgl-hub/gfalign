#include <stdlib.h>
#include <map>
#include <cstdio>

#include "validate.h"

int main(void) {
    std::cout << "WARNING: only run this program if the program is in a working state" << std::endl;
    std::cout << "WARNING: previous validate files will be deleted" << std::endl;
    std::cout << "continue? (Y/N) ";
    std::string input;
    std::cin >> input;
    if(input != "Y" && input != "y") {
        std::cout << "validate generation cancelled" << std::endl;
        std::exit(0);
    }
    std::cout << "deleting old validate files..." << std::endl;

    for(auto &file : list_dir("validateFiles")) {
        if(getFileExt(file) != "tst") continue; // dont delete README
        file = "validateFiles/"+file;
        if(remove(file.c_str()) != 0) {
            std::cerr << "error deleting <" << file << ">" << std::endl;
            return -1;
        }
    }

    std::cout << "generating new validate files..." << std::endl;
    std::vector<std::pair<std::set<std::string>, std::vector<std::string>>> file_args;

    // test eval
    const std::map<std::set<std::string>, std::vector<std::string>> ext_args = {
    //  {{set of test file extensions}, {list of command line args to run with}}
    };
    const std::set<std::string> excludeExt {};
    
    file_args = {
        {{"-f testFiles/random1.gfa"}, {"-g testFiles/random1.gaf", "-g testFiles/random1.gaf --graph-statistics", "-g testFiles/random1.gaf --sort-alignment"}},
        {{"-f testFiles/random2.gfa"}, {"-g testFiles/random2.gaf", "-g testFiles/random2.gaf --graph-statistics", "-g testFiles/random2.gaf --sort-alignment"}}
    //  {{set of test inputs}, {list of command line args to run with}}
    };
    const std::set<std::string> excludeFile {};

    for(const auto &pair : file_args) {
        for(const std::string &input : pair.first) {
            for(const std::string &args : pair.second) {
                genTest("gfalign", "eval", input, args);
            }
        }
    }
    
    // test subgraph
    file_args = {
    //  {{set of test inputs}, {list of command line args to run with}}
    };
    for(const auto &pair : file_args) {
        for(const std::string &input : pair.first) {
            for(const std::string &args : pair.second) {
                genTest("gfalign", "subgraph", input, args);
            }
        }
    }
	
    // test search
    file_args = {
        {{"-f testFiles/random3.gfa"}, {"-n testFiles/random3.search_nodelist.tsv -s 1 -d 4"}}
    //  {{set of test inputs}, {list of command line args to run with}}
    };
    for(const auto &pair : file_args) {
        for(const std::string &input : pair.first) {
            for(const std::string &args : pair.second) {
                genTest("gfalign", "search", input, args);
            }
        }
    }
	
	// test alignment filtering
	file_args = {
		{{"-g testFiles/random3.gaf"}, {"-n testFiles/random3.filter_nodelist.ls -o gaf"}}
	//  {{set of test inputs}, {list of command line args to run with}}
	};
	for(const auto &pair : file_args) {
		for(const std::string &input : pair.first) {
			for(const std::string &args : pair.second) {
				genTest("gfalign", "filter", input, args);
			}
		}
	}
    std::exit(EXIT_SUCCESS);
}
