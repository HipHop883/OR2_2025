#ifndef UTILS_H
#define UTILS_H

#include "tsp.h"

double second();
int compar(const void *a, const void *b);
void set_seed(int seed);
double rand01();
int runs(int argc, char **argv);
void update_perf_csv(const instance *inst, double *run_costs, int num_runs);
void copy_sol(const solution *src, solution *dst);
int is_exact_method(const char *method);

#endif // UTILS_H
