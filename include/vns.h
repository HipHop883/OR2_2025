#ifndef VNS_H_
#define VNS_H_

#include "tsp.h"

int tsp_solve_vns(instance *tsp, solution *sol);

static void generate_3opt_positions(instance *tsp, int *positions);

static int compar(const void *a, const void *b);

int tsp_3opt_solution(instance *tsp, solution *sol);

static void tsp_3opt_swap(instance *tsp, solution *current_sol, int *temp_tour, int *positions);

#endif // VNS_H_
