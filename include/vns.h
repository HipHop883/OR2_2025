#ifndef VNS_H_
#define VNS_H_

#include <tsp.h>

int tsp_solve_vns(instance *tsp, int *output_solution, double *output_cost);

#endif // VNS_H_