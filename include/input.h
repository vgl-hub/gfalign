#ifndef INPUT_H
#define INPUT_H

struct UserInputGfalign : UserInput {
    
    int cmd_flag = 0,
    terminalAlignments_flag = 0,
    sortAlignment_flag = 0;
    std::string nodeList;
    uint32_t dijkstraSteps = 1000;
};

class Input {
    
    UserInputGfalign userInput;
    
    //intermediates
    std::string h;
    
    // stream read variable definition
    std::string firstLine;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
    
public:
    
    void load(UserInputGfalign userInput);
    
    void read(InSequences& inSequence);
    
    void read(InAlignments& inAlignments);
    
};

#endif /* INPUT_H */

