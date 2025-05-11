#ifndef TSP_CPLEX_H_
#define TSP_CPLEX_H_

#include "tsp.h"

int apply_cplex_beneders(instance *inst, solution *sol);
int apply_cplex_branchcut(instance *inst, solution *sol);

#endif /* TSP_CPLEX_H_ */
