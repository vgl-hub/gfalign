#ifndef EVAL_H
#define EVAL_H

void evalGFA(InSequences& InSequences, InAlignments& InAlignments);

void evalPath(InSequences& inSequences, InAlignments& InAlignments, std::string pathStr);

void dijkstra(InSequences& inSequences, InAlignments& InAlignments, std::string nodeFile, std::string source, std::string destination, uint32_t maxSteps, uint32_t minNodes);

#endif /* EVAL_H */
