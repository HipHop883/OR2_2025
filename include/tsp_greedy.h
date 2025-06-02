#ifndef TSP_GREEDY_H_
#define TSP_GREEDY_H_

#include "tsp.h"

/**
 * High-level interface to run multi-start greedy TSP search.
 */
int apply_greedy_search(const instance *inst, solution *sol);

/**
 * Lower-level utility to run nearest neighbor from a fixed start node.
 * Exposed for flexibility and advanced composition use cases.
 */
int apply_nearest_neighbor(const instance *inst, solution *sol, int start_node);

#endif
