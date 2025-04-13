#ifndef TSP_CPLEX_H_
#define TSP_CPLEX_H_

// #include <cplex.h>
#include "tsp.h"

int xpos(int i, int j, const instance *inst);
// int build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
int build_solution(const double *xstar, instance *inst, solution *sol);

// int add_sec(CPXENVptr env, CPXLPptr lp, int *comp, int ncomp, instance *inst);
int build_components(const double *xstar, int *comp, instance *inst);
int apply_cplex_beneders(instance *inst, solution *sol);

#endif /* TSP_CPLEX_H_ */
