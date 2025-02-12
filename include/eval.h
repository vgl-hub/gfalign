#ifndef EVAL_H
#define EVAL_H

void evalGFA(InSequences& InSequences, InAlignments& InAlignments);

void dijkstra(InSequences& inSequences, InAlignments& InAlignments, std::string nodeFile, std::string source, std::string destination, uint32_t maxSteps, int32_t minNodes);

#endif /* EVAL_H */
