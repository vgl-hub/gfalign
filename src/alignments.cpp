#include <stdlib.h>
#include <string>
#include <vector>
#include <mutex>
#include <string.h>

#include <iostream>
#include <fstream>

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "log.h"
#include "global.h"
#include "uid-generator.h"
#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "sak.h"
#include "stream-obj.h"
#include "input-agp.h"
#include "input-filters.h"
#include "input-gfa.h"

#include "alignments.h"

InAlignment::InAlignment(std::vector<std::string> cols, std::vector<Tag> inTags, unsigned int pos) {
    
    this->qName = cols[0];
    this->qLen = stoi(cols[1]);
    this->qStart = stoi(cols[2]);
    this->qEnd = stoi(cols[3]);
    this->strand = cols[4][0];
    this->path = cols[5];
    this->pLen = stoi(cols[6]);
    this->pStart = stoi(cols[7]);
    this->pEnd = stoi(cols[8]);
    this->matches = stoi(cols[9]);
    this->blockLen = stoi(cols[10]);
    this->mapq = stoi(cols[11]);
    this->inTags = inTags;
    
    this->pos = pos;

}

std::string InAlignment::print() {
    
    std::string alignment =
    qName + "\t" +
    std::to_string(qLen) + "\t" +
    std::to_string(qStart) + "\t" +
    std::to_string(qEnd) + "\t" +
    std::string(&strand, 1) + "\t" +
    path + "\t" +
    std::to_string(pLen) + "\t" +
    std::to_string(pStart) + "\t" +
    std::to_string(pEnd) + "\t" +
    std::to_string(matches) + "\t" +
    std::to_string(blockLen) + "\t" +
    std::to_string(mapq);
    
    for (Tag tag : inTags) {
    
        alignment += std::string("\t") + tag.label + std::string(":") + tag.type + std::string(":") + tag.content;
        
    }
    
    alignment += "\n";
    
    return alignment;

}

InAlignments::~InAlignments() {

    for (InAlignment* p : inAlignments)
        delete p;
}

void InAlignments::load(std::shared_ptr<std::istream> stream, int terminalAlignments_flag) {

    this->terminalAlignments_flag = terminalAlignments_flag;
    unsigned int batchSize = 10000;
    StreamObj streamObj;
    std::string* alignment = new std::string;
    Alignments* alignmentBatch = new Alignments;
    std::vector<std::string> arguments;

    if (stream) {

        while (getline(*stream, *alignment)) {

            arguments = readDelimited(*alignment, "\t");

            alignmentBatch->alignments.push_back(alignment);
            pos++;
            
            alignment = new std::string;

            if (pos % batchSize == 0) {
                
                if (stream->eof())
                    break;
                
                alignmentBatch->batchN = pos/batchSize;
                
                lg.verbose("Processing batch N: " + std::to_string(alignmentBatch->batchN));

                appendAlignments(alignmentBatch);
                
                alignmentBatch = new Alignments;

            }

        }
        
        alignmentBatch->batchN = pos/batchSize+1;
        
        lg.verbose("Processing batch N: " + std::to_string(alignmentBatch->batchN));

        appendAlignments(alignmentBatch);

    }

}

void InAlignments::appendAlignments(Alignments* alignmentBatch) { // read a collection of alignments
    
    threadPool.queueJob([=]{ return traverseInAlignments(alignmentBatch); });
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    lck.unlock();
    
}

bool InAlignments::traverseInAlignments(Alignments* alignmentBatch) { // traverse the read

    Log threadLog;
    
    threadLog.setId(alignmentBatch->batchN);
    
    std::vector<InAlignment*> inAlignmentsBatch;
    
    AlignmentStats tmpStats;
    
    unsigned int readN = 0;
    
    for (std::string* alignment : alignmentBatch->alignments) {
        
        inAlignmentsBatch.push_back(traverseInAlignment(&threadLog, alignment, alignmentBatch->batchN+readN++, &tmpStats));
        
    }
    
    delete alignmentBatch;
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    updateStats(&tmpStats);
    
    inAlignments.insert(std::end(inAlignments), std::begin(inAlignmentsBatch), std::end(inAlignmentsBatch));
    
    logs.push_back(threadLog);
    
    lck.unlock();
    
    return true;
    
}

InAlignment* InAlignments::traverseInAlignment(Log* threadLog, std::string* alignment, unsigned int pos, AlignmentStats* tmpStats) { // traverse a single read
    
    std::vector<std::string> cols = readDelimited(*alignment, "\t");
    
    std::vector<Tag> inTags;
    
    std::vector<std::string> tagValues;
    
    Tag tag;
    
    for (unsigned int i = 12; i < cols.size(); i++) {
        
        tagValues = readDelimited(cols[i], ":");
        
        tag.label[0] = tagValues[0][0];
        tag.label[1] = tagValues[0][1];
        tag.type = tagValues[1][0];
        tag.content = tagValues[2];
    
        inTags.push_back(tag);
    
    }

    InAlignment* inAlignment = new InAlignment(cols, inTags, pos);
    
    delete alignment;
    
    tmpStats->add(inAlignment);
    
    threadLog->add("Individual alignment read: " + cols[0]);
    
    return inAlignment;
    
}

void AlignmentStats::add(InAlignment* alignment){
    
    tmpQLen += alignment->qLen;
    tmpAlgSeq += alignment->qEnd - alignment->qStart;
    alignment->strand == '+' ? ++plus : ++minus;
    tmpPLen += alignment->pLen;
    tmpMatches += alignment->matches;
    tmpBlockLen += alignment->blockLen;
    tmpMapq += alignment->mapq;
    
}

void InAlignments::updateStats(AlignmentStats* tmpStats) {
    
    totQLen += tmpStats->tmpQLen;
    totAlgSeq += tmpStats->tmpAlgSeq;
    totPlus += tmpStats->plus;
    totMinus += tmpStats->minus;
    totPLen += tmpStats->tmpPLen;
    totMatches += tmpStats->tmpMatches;
    totBlockLen += tmpStats->tmpBlockLen;
    totMapq += tmpStats->tmpMapq;
    
}

void InAlignments::printStats() {
    
    if (!tabular_flag) {
    
        std::cout<<output("+++Alignment summary+++")<<"\n";
    
    }

    std::cout<<output("# alignments")<<inAlignments.size()<<"\n";
    std::cout<<output("Average read length")<<gfa_round(computeAvg(totQLen))<<"\n";
    std::cout<<output("Average aligned sequence")<<gfa_round(computeAvg(totAlgSeq))<<"\n";
    std::cout<<output("Alignment orientation (+/-)")<<totPlus<<"("<<gfa_round((double) totPlus/(totPlus+totMinus)*100)<<"%):"<<totMinus<<"("<<gfa_round((double) totMinus/(totPlus+totMinus)*100)<<"%)"<<"\n";
    std::cout<<output("Average path length")<<gfa_round(computeAvg(totPLen))<<"\n";
    std::cout<<output("Average alignment quality")<<gfa_round(computeAvg(totMapq))<<"\n";
    std::cout<<output("Average matches #")<<gfa_round(computeAvg(totMatches))<<"\n";
    std::cout<<output("Average block length")<<gfa_round(computeAvg(totBlockLen))<<"\n";
    std::cout<<output("Primary alignments")<<primaryAlignments<<"\n";
    std::cout<<output("Secondary alignments")<<secondaryAlignments<<"\n";
    std::cout<<output("Supplementary alignments")<<supplementaryAlignments<<"\n";
    std::cout<<output("Terminal supplementary alignments")<<terminalSupplementaryAlignments<<"\n";

}

double InAlignments::computeAvg(long long unsigned int value) {
    
    return (double) value/inAlignments.size();
    
}

void InAlignments::sortAlignmentsByNameAscending(){
    
    sort(inAlignments.begin(), inAlignments.end(), [](InAlignment* one, InAlignment* two){return one->qName < two->qName;});
    
}

void InAlignments::outAlignments(){
    
    for(InAlignment* alignment : inAlignments) {
        
        std::cout<<alignment->print();
        
    }
    
}

void InAlignments::markDuplicates(){
    
    std::string prevQname;
    
    std::vector<InAlignment*> alignments;
    
    for(auto alignment = inAlignments.begin(); alignment != inAlignments.end(); ++alignment) {
        
        alignments.push_back(*alignment);
        
        if ((*alignment)->qName == prevQname) {
            
            ++secondaryAlignments;
            
            if(std::next(alignment) == inAlignments.end() || (*std::next(alignment))->qName != (*alignment)->qName){
                
                countSupplementary(alignments);
                
                alignments.clear();
                
            }
            
        }else{
            
            ++primaryAlignments;
            
            prevQname = (*alignment)->qName;
            
        }
        
    }
    
}

void InAlignments::countSupplementary(std::vector<InAlignment*> alignments){
    
    unsigned int pos = 0, count = 0;
    
    sort(alignments.begin(), alignments.end(), [](InAlignment* one, InAlignment* two){return one->qStart < two->qStart;}); // first sort the alignments by their qStart
    
    for (InAlignment* alignment : alignments) {
    
        if(pos != 0 && alignment->qStart > pos) { // if this is not the first alignment and we are aligning a downstream portion of the read
            
            ++supplementaryAlignments;
            
            ++count;
            
        }
        
        pos = alignment->qEnd; // ensure alignments do not overlap
        
    }
    
    if (alignments.size() == 2 && count == 1) { // we are only looking at unambigous supplementary alignments
        
        if (alignments[0]->pEnd >= alignments[0]->pLen - 500 && alignments[1]->pStart <= 500) { // the end of the leftmost alignment ends at the path end and the start of the rightmost alignment is at the beginning of the path
            
            ++terminalSupplementaryAlignments;
            
            if (terminalAlignments_flag)
                std::cout<<alignments[0]->print()<<alignments[1]->print();
        }
        
    }
    
}

void InAlignments::buildEdgeGraph(phmap::flat_hash_map<std::string, unsigned int>* headersToIds, phmap::flat_hash_map<unsigned int, std::string>* idsToHeaders, unsigned int uId) { // graph constructor
    
    lg.verbose("Started edge graph construction from alignment");
    
    adjEdgeList.clear();
    
    adjEdgeList.resize(uId); // resize the adjaciency list to hold all nodes
    
    for (InAlignment* alignment : inAlignments) // search candidate edges in the alignment
    {
        
        std::vector<InEdge> edges = GAFpathToEdges(alignment->path, headersToIds);
        
        for (auto &edge: edges) // add edges to the graph
        {
            
            Edge fwEdge {edge.getsId1Or(), edge.getsId2(), edge.getsId2Or()};
            
            auto it = find(adjEdgeList.at(edge.getsId1()).begin(), adjEdgeList.at(edge.getsId1()).end(), fwEdge);
            
            if (it == adjEdgeList.at(edge.getsId1()).end()) {
            
                lg.verbose("Adding edge: " + (*idsToHeaders)[edge.getsId1()] + "(" + std::to_string(edge.getsId1()) + ") " + edge.getsId1Or() + " " + (*idsToHeaders)[edge.getsId2()] + "(" + std::to_string(edge.getsId2()) + ") " + edge.getsId2Or());

                adjEdgeList.at(edge.getsId1()).push_back({edge.getsId1Or(), edge.getsId2(), edge.getsId2Or(), 1}); // insert at edge start gap destination and orientations

                Edge rvEdge {edge.getsId2Or() == '+' ? '-' : '+', edge.getsId1(), edge.getsId1Or() == '+' ? '-' : '+', 1};

                if (find(adjEdgeList.at(edge.getsId2()).begin(), adjEdgeList.at(edge.getsId2()).end(), rvEdge) == adjEdgeList.at(edge.getsId2()).end()) // add backward edge only if is not already present
                    adjEdgeList.at(edge.getsId2()).push_back(rvEdge); // assembly are bidirected by definition
            
            }else{
                
                lg.verbose("Edge already present, increasing weight: " + (*idsToHeaders)[edge.getsId1()] + "(" + std::to_string(edge.getsId1()) + ") " + edge.getsId1Or() + " " + (*idsToHeaders)[edge.getsId2()] + "(" + std::to_string(edge.getsId2()) + ") " + edge.getsId2Or());
                
                it->weight++;
                
                Edge rvEdge {edge.getsId2Or() == '+' ? '-' : '+', edge.getsId1(), edge.getsId1Or() == '+' ? '-' : '+', 1};
                
                auto it2 = find(adjEdgeList.at(edge.getsId2()).begin(), adjEdgeList.at(edge.getsId2()).end(), rvEdge);
                
                it2->weight++;
                
            }
            
        }
    }
    
    lg.verbose("Graph built");
    
}

std::vector<std::vector<Edge>> InAlignments::getEdgeGraph() {
    
    return adjEdgeList;
    
}

std::vector<InEdge> GAFpathToEdges(std::string path, phmap::flat_hash_map<std::string, unsigned int>* headersToIds) {

    std::vector<InEdge> edges;
    size_t pos = 0;
    std::string sId1, sId2;
    char sId1Or = '+', sId2Or = '+';
    unsigned int counter = 0;

    while (path.size() != 0) {
        
        if(path[pos] == '>' || path[pos] == '<' || pos == path.size()) {
            
            if (pos == 0) {
                pos++;
                continue;
            }
            
            counter++;
            
            (counter == 1 ? sId1Or : sId2Or) = (path[0] == '>' ? '+' : '-');
            
            (counter == 1 ? sId1 : sId2) = path.substr(1, pos - 1);
            
            path.erase(0, pos);
            
            if (counter == 2) {
                
                InEdge edge;
                
                edge.newEdge(0, (*headersToIds)[sId1], (*headersToIds)[sId2], sId1Or, sId2Or);
                
                edges.push_back(edge);
                
                sId1Or = sId2Or;
                sId1 = sId2;
                
                counter = 1;
                
            }
            
            pos = 0;
            
        }else{
    
            pos++;
        
        }
            
    }
        
    return edges;
    
}
