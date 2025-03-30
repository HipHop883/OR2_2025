#ifndef TSP_HEURISTICS_H_
#define TSP_HEURISTICS_H_

#include "tsp.h"
#include "tsp_greedy.h"

// data structure for the tabu list
typedef struct
{
    int **tabu_list;
    int tenure;
    int size;
} tabuList;

int apply_heuristic_vns(instance *tsp, solution *sol);
int apply_heuristic_tabu(instance *inst, solution *sol);

#endif // TSP_HEURISTICS_H_
