#ifndef ALIGNMENTS_H
#define ALIGNMENTS_H

class InAlignment{
    
    std::string qName;
    unsigned int qLen;
    unsigned int qStart;
    unsigned int qEnd;
    char strand;
    std::string path;
    unsigned int pLen;
    unsigned int pStart;
    unsigned int pEnd;
    unsigned int matches;
    unsigned int blockLen;
    unsigned int mapq;
    unsigned int pos;

public:
    
    InAlignment(std::vector<std::string> cols, unsigned int pos);
    std::string print();
    
};

class InAlignments{
    
    std::vector<Log> logs;
    
    UserInput userInput;
    std::vector<InAlignment*> inAlignments;
    
    unsigned int pos = 0; // to keep track of the original order
    
public:
    
    ~InAlignments();
    
    void load(UserInput userInput);
    
    void traverseInAlignments(Alignments* sequence);
    
    InAlignment* traverseInAlignment(Log* threadLog, std::string* alignment, unsigned int pos);
    
    void appendAlignments(Alignments* alignmentBatch);
    
    unsigned long long int getTotReadLen();
    
    double computeAvgReadLen();
    
    unsigned long long int getReadN50();
    
    void report();
    
    void evalNstars();
    
};

#endif /* ALIGNMENTS_H */

