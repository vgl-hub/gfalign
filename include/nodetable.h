#ifndef NODETABLE_H
#define NODETABLE_H

struct Record {
	uint32_t uId, count;
	//char type;
};

struct NodeTable {
	
	phmap::flat_hash_map<std::string,Record> records;
	uint32_t nodeCount = 0;
	
	NodeTable() {}
	
	NodeTable(std::string nodeFile, phmap::flat_hash_map<std::string,uint32_t> &lookupTable) {
		
		std::string line;
		std::ifstream file(nodeFile);
		while (std::getline(file, line)) {
			
			std::vector<std::string> lineVec = readDelimited(line, "\t");
			
			auto got = lookupTable.find(lineVec.at(0)); // find uId for the node
			uint32_t count = 1;
			if (lineVec.size() > 1) { // the file contains a count column
				count = std::stoi(lineVec.at(1));
				if (count < 1) // if count is 0 don't introduce the record
					continue;
			}
			nodeCount += count; // keep track of the total number of nodes for a Hamiltonian path
			if (got != lookupTable.end()) {
				Record record{got->second,count};
				records.insert(std::make_pair(lineVec.at(0),record));
				lg.verbose("Added record: " + lineVec.at(0) + " " + std::to_string(record.count));
			}else{
				fprintf(stderr, "Error: node not in graph (pIUd: %s)\n", lineVec.at(0).c_str());
				exit(EXIT_FAILURE);
			}
		}
		file.close();
		lg.verbose("Node table read");
	}
	
	Record operator [](std::string node) const {
		return records.at(node);
	}
	
	void add(std::string node, Record record) {
		if (record.count < 1) // only insert record if the count is > 0
			return;
		records.insert(std::make_pair(node,record));
		nodeCount += record.count;
	}
	
	bool checkHamiltonian(phmap::flat_hash_map<uint32_t,uint32_t> &pathNodes, uint32_t pathNodesCount) {
		if (pathNodesCount + 2 != nodeCount)
			return false;
		for (auto& it: records) {
			auto found = pathNodes.find(it.second.uId);
			if (found == pathNodes.end() || found->second != it.second.count) {
				lg.verbose("This is not a Hamiltonian path.");
				return false;
			}
		}
		return true;
	}
};

#endif // #ifndef NODETABLE_H
