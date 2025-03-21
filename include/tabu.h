#ifndef TABU_H_
#define TABU_H_

#include <tsp.h>

int tabu_method(instance *inst, solution *sol);

static int two_opt_random(instance *tsp, solution *sol, int **tabu_list);

#endif //TABU_H_