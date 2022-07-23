#ifndef INPUT_H
#define INPUT_H

class Input {
    
    UserInput userInput;
    
    //intermediates
    std::string h;
    char* c;
    
    // stream read variable definition
    std::string firstLine;
    bool stopStream = false;
    unsigned int pos = 0; // to keep track of the original order
    
    StreamObj streamObj;
    
    std::string* alignment;
    
    std::shared_ptr<std::istream> stream;
    
public:
    
    void load(UserInput userInput);
    
    void read(InSequences& inSequence);
    
    void read(InAlignments& inAlignments);
    
};

#endif /* INPUT_H */

