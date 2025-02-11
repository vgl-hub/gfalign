#ifndef ALIGNMENTS_H
#define ALIGNMENTS_H

struct Alignments { // collection of sequences
    
    std::vector<std::string*> alignments;
    unsigned int batchN;
    
};

struct Step { // step in an alignment path
	int32_t id;
	char orientation;
	
	inline bool operator==(const Step& step) const {
		return this->id == step.id && this->orientation == step.orientation;
	}
	bool operator !=(const Step& step) const {
		return this->id != step.id && this->orientation != step.orientation;
	}
};

struct Path { // graph alignment path
	std::vector<Step> path;
	NodeTable nodeTable;
	
	Path() {}
	
	Path(NodeTable nodeTable) : nodeTable(nodeTable) {}
	
	uint32_t size() const {
		return path.size();
	}
	
	void push_back(int32_t id, char orientation) {
		path.push_back(Step{id, orientation});
	}
	
	const Step& operator[](size_t index) const {
		if (index >= size()) {
			throw std::out_of_range("Index out of bounds");
		}
		return path[index];
	}
	
	const Step& at(size_t index) const {
		if (index >= size()) {
			throw std::out_of_range("Index out of bounds");
		}
		return path.at(index);
	}
	
	std::unordered_set<uint32_t> pathToSet() const {
		std::unordered_set<uint32_t> uIdSet;
		for (Step step : path)
			uIdSet.insert(step.id);
		return uIdSet;
	}
	
	void reverse() {
		std::reverse(path.begin(), path.end());
	}
	
	void print() const {
		for (uint32_t i = 0; i<path.size(); ++i) {
			std::cout<<path.at(i).id<<path.at(i).orientation;
			if (i+1<path.size())
				std::cout<<',';
		}
		std::cout<<std::endl;
	}
};

struct PairwisePathAlignment{
	
	Path A, B;
	
	PairwisePathAlignment(Path A, Path B) : A(A), B(B){}
	
	void print(phmap::flat_hash_map<unsigned int, std::string> &idsToHeaders) const {
		
		for (uint32_t i = 0; i<A.size(); ++i) {
			if (A.at(i).id == -1)
				std::cout<<std::string(idsToHeaders[B.at(i).id].size()+1, '-')<<",";
			else
				std::cout<<idsToHeaders[A.at(i).id]<<A.at(i).orientation<<",";
		}
		std::cout<<std::endl;
		for (uint32_t i = 0; i<B.size(); ++i) {
			if (B.at(i).id == -1)
				std::cout<<std::string(idsToHeaders[A.at(i).id].size()+1, '-')<<",";
			else
				std::cout<<idsToHeaders[B.at(i).id]<<B.at(i).orientation<<",";
		}
		std::cout<<std::endl;
	}
};

class InAlignment{ // full graph alignment record (GAF)
    
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
	
	Path GAFpathToPath(phmap::flat_hash_map<std::string, unsigned int> &headersToIds);
    
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

class InAlignments{ // collection of alignments
    
    int terminalAlignments_flag = 0;
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
                            totBlockLen = 0,
                            primaryAlignments = 0,
                            secondaryAlignments = 0,
                            supplementaryAlignments = 0,
                            terminalSupplementaryAlignments = 0;
    
    std::vector<std::vector<Edge>> adjEdgeList;
    
public:
    void buildEdgeGraph(phmap::flat_hash_map<std::string, unsigned int>* headersToIds, phmap::flat_hash_map<unsigned int, std::string>* idsToHeaders, unsigned int uId);
    
public:
    
    ~InAlignments();
    
    void load(std::string file, int terminalAlignments_flag);
    
    bool traverseInAlignments(Alignments* sequence);
    
    InAlignment* traverseInAlignment(Log* threadLog, std::string* alignment, unsigned int pos, AlignmentStats* tmpStats);
    
    void appendAlignments(Alignments* alignmentBatch);
    
    void printStats();
    
    unsigned long long int getTotAlignments();
    
    void updateStats(AlignmentStats* tmpStats);
    
    double computeAvg(long long unsigned int value);
    
    void sortAlignmentsByNameAscending();
    
    void markDuplicates();
    
    void outAlignments();
    
    std::vector<std::vector<Edge>> getEdgeGraph();
    
    void countSupplementary(std::vector<InAlignment*> alignments);
    
    void sortAlignmentsByQStart(std::vector<InAlignment*>* alignments);
	
	std::vector<InAlignment*> getAlignments() const;
	
	std::vector<Path> getPaths(phmap::flat_hash_map<std::string, unsigned int> &headersToIds);
    
};

std::vector<InEdge> GAFpathToEdges(std::string path, phmap::flat_hash_map<std::string, unsigned int>* headersToIds);

#define MAX_N 1001
PairwisePathAlignment alignPaths(int8_t match_score, int8_t mismatch_score, int8_t gap_score, Path A, Path B, int dp[MAX_N][MAX_N]);

#endif /* ALIGNMENTS_H */

