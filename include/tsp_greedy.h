#ifndef TSP_GREEDY_H_
#define TSP_GREEDY_H_

#include "tsp.h"

int apply_greedy_search(const instance *inst, solution *sol);
int apply_nearest_neighbor(const instance *inst, solution *sol, int start_node);

#endif /* TSP_GREEDY_H_ */
