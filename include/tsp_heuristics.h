#ifndef TSP_HEURISTICS_H_
#define TSP_HEURISTICS_H_

#include "tsp.h"

#define NO_IMPROVE_LIMIT 50

// data structure for the tabu list
typedef struct
{
    int **tabu_list;
    int tenure;
    int size;

    // For dynamic tenure
    int min_tenure;
    int max_tenure;
    int iterations_without_improvement;
} tabuList;

int apply_heuristic_vns(instance *tsp, solution *sol);
int apply_heuristic_tabu(instance *inst, solution *sol);

#endif // TSP_HEURISTICS_H_
