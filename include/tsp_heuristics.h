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

int tsp_solve_vns(instance *tsp, solution *sol);
int tsp_solve_tabu(instance *inst, solution *sol);

static int best_2opt_not_tabu(instance *tsp, solution *sol, tabuList *tabu);
int is_tabu_move(tabuList *tabu, int i, int j);
void add_tabu_move(tabuList *tabu, int i, int j);
void free_tabu(tabuList *tabu);
tabuList *init_tabu_list(int nnodes);

#endif // TSP_HEURISTICS_H_
