#ifndef TSP_GREEDY_H_
#define TSP_GREEDY_H

#include "tsp.h"

int nearest_neighbor(const instance *inst, solution *sol, int start_node);
int solve_greedy(const instance *inst, solution *sol);

#endif
