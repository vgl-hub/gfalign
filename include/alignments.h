#ifndef ALIGNMENTS_H
#define ALIGNMENTS_H

struct Alignments { // collection of alignments
    
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
		return this->id != step.id || this->orientation != step.orientation;
	}
};

struct Path { // graph alignment path
	std::string qName;
	std::vector<Step> path;
	NodeTable nodeTable;
	
	Path() {}
	
	Path(std::string qName) : qName(qName) {}
	
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
		return path.at(index);
	}
	
	phmap::flat_hash_map<uint32_t,uint32_t> pathToMap() const {
		phmap::flat_hash_map<uint32_t,uint32_t> uIdMap;
		for (Step step : path)
			++uIdMap[step.id];
		return uIdMap;
	}
	
	void reverse() {
		std::reverse(path.begin(), path.end());
	}
	
	Path reverseComplement() {
		Path rc = *this;
		std::reverse(rc.path.begin(), rc.path.end());
		for (Step &step : rc.path)
			step.orientation = (step.orientation == '+') ? '-' : '+';
		return rc;
	}
	
	void print(InSequences &inSequences) const {
		for (uint32_t i = 0; i<this->path.size(); ++i) {
			std::cout<<inSequences.findSegmentBySUId(this->path.at(i).id).getSeqHeader()<<this->path.at(i).orientation;
			if (i+1<this->path.size())
				std::cout<<',';
		}
	}
};

struct PairwisePathAlignment{
	
	Path A, B;
	int32_t alignmentScore = 0;
	
	PairwisePathAlignment(Path A, Path B, int32_t alignmentScore) : A(A), B(B), alignmentScore(alignmentScore) {}
	
	std::string getAlignment(phmap::flat_hash_map<unsigned int, std::string> &idsToHeaders, bool doNotReturnRef = false) const {
		
		std::string aln;
		
		if (!doNotReturnRef) {
			for (uint32_t i = 0; i<A.size(); ++i) {
				if (A.at(i).id == -1)
					aln += std::string(idsToHeaders[B.at(i).id].size()+1, '-') + ",";
				else if (A.at(i) != B.at(i))
					aln += idsToHeaders[A.at(i).id] + A.at(i).orientation + ",";
				else
					aln += std::string(idsToHeaders[A.at(i).id].size()+1, '.') + ",";
			}
			aln += '\n';
		}
		for (uint32_t i = 0; i<B.size(); ++i) {
			if (B.at(i).id == -1)
				aln += std::string(idsToHeaders[A.at(i).id].size()+1, '-') + ",";
			else if (A.at(i) != B.at(i))
				aln += idsToHeaders[B.at(i).id] + B.at(i).orientation + ",";
			else
				aln += std::string(idsToHeaders[B.at(i).id].size()+1, '.') + ",";
		}
		return aln;
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
	
	bool isContained(phmap::flat_hash_set<std::string> &headers);
	
	uint32_t pathNodesCount();
    
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
    
    void outputAlignments(std::string file);
    
    std::vector<std::vector<Edge>> getEdgeGraph();
    
    void countSupplementary(std::vector<InAlignment*> alignments);
    
    void sortAlignmentsByQStart(std::vector<InAlignment*>* alignments);
	
	std::vector<InAlignment*> getAlignments() const;
	
	std::vector<Path> getPaths(phmap::flat_hash_map<std::string, unsigned int> &headersToIds);
	
	void filterAlignmentByNodelist(std::vector<std::string> nodelist, int32_t minNodes);
    
};

std::vector<InEdge> GAFpathToEdges(std::string path, phmap::flat_hash_map<std::string, unsigned int>* headersToIds);

#define MAX_N 1001
PairwisePathAlignment alignPaths(int8_t match_score, int8_t mismatch_score, int8_t gap_score, Path A, Path B, int dp[MAX_N][MAX_N]);

#endif /* ALIGNMENTS_H */

