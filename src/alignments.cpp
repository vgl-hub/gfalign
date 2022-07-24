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

#include "bed.h"
#include "struct.h"
#include "gfa-lines.h"
#include "stream-obj.h"

#include "functions.h" // global functions

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
    std::to_string(strand) + "\t" +
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

InAlignments::~InAlignments()
{

//    for (InSegment* p : inAlignments)
//        delete p;

}

void InAlignments::load(UserInput userInput) {

    unsigned int batchSize = 10000;
    
    StreamObj streamObj;
    
    std::string* alignment = new std::string;
    
    std::shared_ptr<std::istream> stream;

    stream = streamObj.openStream(userInput, 'g');

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

void InAlignments::traverseInAlignments(Alignments* alignmentBatch) { // traverse the read

    Log threadLog;
    
    threadLog.setId(alignmentBatch->batchN);
    
    std::vector<InAlignment*> inAlignmentsBatch;
    
    unsigned int readN = 0;
    
    for (std::string* alignment : alignmentBatch->alignments) {
        
        inAlignmentsBatch.push_back(traverseInAlignment(&threadLog, alignment, alignmentBatch->batchN+readN++));
        
    }
    
    delete alignmentBatch;
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    inAlignments.insert(std::end(inAlignments), std::begin(inAlignmentsBatch), std::end(inAlignmentsBatch));
    
    logs.push_back(threadLog);
    
    lck.unlock();
    
}

InAlignment* InAlignments::traverseInAlignment(Log* threadLog, std::string* alignment, unsigned int pos) { // traverse a single read
    
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
    
    threadLog->verbose("Individual alignment read: " + cols[0]);
    
    return inAlignment;
    
}

void InAlignments::printStats() {

    std::cout<<output("# alignments")<<getTotAlignments()<<"\n";
    
    computeStats();
    
    std::cout<<output("Average alignment quality")<<getAvgQual()<<"\n";

}

void InAlignments::computeStats() {
    
    unsigned long long int totAlignments = inAlignments.size();
    
    for(InAlignment* alignment : inAlignments) {
        
        avgQual += alignment->mapq;
        
    }
    
    avgQual = avgQual / totAlignments;
    
}

unsigned long long int InAlignments::getTotAlignments() {
    
    return inAlignments.size();
    
}

unsigned int InAlignments::getAvgQual() {
    
    return avgQual;
    
}