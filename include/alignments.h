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
    std::vector<Tag> inTags;
    
    unsigned int pos;

public:
    
    InAlignment(std::vector<std::string> cols, std::vector<Tag> inTags, unsigned int pos);
    std::string print();
    
    friend class InAlignments;
    
};

class InAlignments{
    
    std::vector<Log> logs;
    
    UserInput userInput;
    std::vector<InAlignment*> inAlignments;
    
    unsigned int pos = 0;
    
    unsigned int avgQual;
    
public:
    
    ~InAlignments();
    
    void load(UserInput userInput);
    
    void traverseInAlignments(Alignments* sequence);
    
    InAlignment* traverseInAlignment(Log* threadLog, std::string* alignment, unsigned int pos);
    
    void appendAlignments(Alignments* alignmentBatch);
    
    void printStats();
    
    unsigned long long int getTotAlignments();
    
    void computeStats();
    
    unsigned int getAvgQual();
    
    
    
};

#endif /* ALIGNMENTS_H */

