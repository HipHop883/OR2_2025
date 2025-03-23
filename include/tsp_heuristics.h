#ifndef TSP_HEURISTICS_H_
#define TSP_HEURISTICS_H_

#include "tsp.h"
#include "tsp_greedy.h"

int tsp_solve_vns(instance *tsp, solution *sol);
int tsp_solve_tabu(instance *inst, solution *sol);

static int two_opt_random(instance *tsp, solution *sol, int **tabu_list);

#endif // TSP_HEURISTICS_H_
