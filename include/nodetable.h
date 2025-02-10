#ifndef NODETABLE_H
#define NODETABLE_H

struct Record {
	uint32_t uId, count;
	//char type;
};

struct NodeTable {
	
	phmap::flat_hash_map<std::string,Record> records;
	
	NodeTable() {}
	
	NodeTable(std::string nodeFile, phmap::flat_hash_map<std::string,uint32_t> &lookupTable) {
		
		std::string line;
		std::ifstream file(nodeFile);
		while (std::getline(file, line)) {
			
			std::vector<std::string> lineVec = readDelimited(line, "\t");
			
			auto got = lookupTable.find(lineVec.at(0));
			uint32_t count = 1;
			if (lineVec.size() > 1)
				count = std::stoi(lineVec.at(1));
			
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
		records.insert(std::make_pair(node,record));
	}
};

#endif // #ifndef NODETABLE_H
