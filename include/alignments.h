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
    
    friend class InAlignments;
    friend class AlignmentStats;
    
};

class AlignmentStats {
    
    long long unsigned int  tmpQLen = 0,
                            tmpAlgSeq = 0, // total aligned sequence
                            plus = 0, minus = 0,
                            tmpPLen = 0,
                            tmpMapq = 0,
                            tmpMatches = 0,
                            tmpBlockLen = 0;

public:
    void add(InAlignment* alignment);
    
    friend class InAlignments;
        
};

class InAlignments{
    
    std::vector<Log> logs;
    
    UserInput userInput;
    std::vector<InAlignment*> inAlignments;
    
    unsigned int pos = 0;
    
    long long unsigned int  totQLen = 0,
                            totAlgSeq = 0,
                            totPlus = 0, totMinus = 0,
                            totPLen = 0,
                            totMapq = 0,
                            totMatches = 0,
                            totBlockLen = 0;
    
    std::vector<std::vector<Edge>> adjEdgeList;
    
public:
    void buildEdgeGraph(phmap::flat_hash_map<std::string, unsigned int>* headersToIds, phmap::flat_hash_map<unsigned int, std::string>* idsToHeaders, unsigned int uId);
    
public:
    
//    ~InAlignments();
    
    void load(UserInput userInput);
    
    bool traverseInAlignments(Alignments* sequence);
    
    InAlignment* traverseInAlignment(Log* threadLog, std::string* alignment, unsigned int pos, AlignmentStats* tmpStats);
    
    void appendAlignments(Alignments* alignmentBatch);
    
    void printStats();
    
    unsigned long long int getTotAlignments();
    
    void updateStats(AlignmentStats* tmpStats);
    
    double computeAvg(long long unsigned int value);
    
    std::vector<std::vector<Edge>> getEdgeGraph();
    
    
};


#endif /* ALIGNMENTS_H */

