#ifndef INPUT_H
#define INPUT_H

struct UserInputGfalign : UserInput {
    
    int32_t cmd_flag = 0,
    alignStats_flag = 0,
    terminalAlignments_flag = 0,
	sortAlignment_flag = 0;
	uint32_t minNodes = 0;
    std::string nodeFile, source, destination;
    uint32_t dijkstraSteps = 100000;
};

class Input {
    
    UserInputGfalign userInput;
	InSequences inSequences;
	InAlignments inAlignments;
    
    //intermediates
    std::string h;
    
    // stream read variable definition
    std::string firstLine;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
    
public:
    
    void loadInput(UserInputGfalign userInput);
    
    void read();
    
    void read(InAlignments& inAlignments);
	
	std::vector<std::string> readNodelist();
    
};

#endif /* INPUT_H */

