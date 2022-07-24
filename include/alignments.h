#ifndef ALIGNMENTS_H
#define ALIGNMENTS_H

struct Alignments { // a collection of sequences
    
    std::vector<std::string*> alignments;
    unsigned int batchN;
    
};

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
    
    unsigned int getMatches();
    unsigned int getBlockLen();
    unsigned int getMapq();
    
    friend class InAlignments;
    
};

struct Stats {
    
    unsigned int tmpMatches = 0;
    unsigned int tmpBlockLen = 0;
    unsigned int tmpMapq = 0;
    
    void add(InAlignment* alignment){
        
        tmpMatches += alignment->getMatches();
        tmpBlockLen += alignment->getBlockLen();
        tmpMapq += alignment->getMapq();
        
    }
        
};

class InAlignments{
    
    std::vector<Log> logs;
    
    UserInput userInput;
    std::vector<InAlignment*> inAlignments;
    
    unsigned int pos = 0;
    
    unsigned int totMatches = 0;
    unsigned int totBlockLen = 0;
    unsigned int totMapq = 0;
    
public:
    
//    ~InAlignments();
    
    void load(UserInput userInput);
    
    void traverseInAlignments(Alignments* sequence);
    
    InAlignment* traverseInAlignment(Log* threadLog, std::string* alignment, unsigned int pos, Stats* tmpStats);
    
    void appendAlignments(Alignments* alignmentBatch);
    
    void printStats();
    
    unsigned long long int getTotAlignments();
    
    void updateStats(Stats* tmpStats);
    
    double getAvgQual();
    
    double getAvgMatches();
    
    double getAvgBlockLen();
    
    
};


#endif /* ALIGNMENTS_H */

