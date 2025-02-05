#ifndef EVAL_H
#define EVAL_H

void evalGFA(InSequences& InSequences, InAlignments& InAlignments);

void dijkstra(InSequences& inSequences, std::vector<std::string> nodeList, std::string source, std::string destination, uint32_t maxSteps);

#endif /* EVAL_H */
