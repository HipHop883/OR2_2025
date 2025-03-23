#ifndef TSP_HEURISTICS_H_
#define TSP_HEURISTICS_H_

#include "tsp.h"
#include "tsp_greedy.h"

#define MIN_TENURE 2
#define MAX_TENURE 15

// data structure for the tabu list
typedef struct{
    int **tabu_list;
    int tenure; 
    int size;
} tabuList;


int tsp_solve_vns(instance *tsp, solution *sol);
int tsp_solve_tabu(instance *inst, solution *sol);

static int two_opt_random(instance *tsp, solution *sol, tabuList *tabu);
int check_tabu_list(tabuList *tabu, int i, int j);
int increase_list(tabuList *tabu);
int decrease_list(tabuList *tabu);


#endif // TSP_HEURISTICS_H_
