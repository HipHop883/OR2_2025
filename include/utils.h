#ifndef UTILS_H
#define UTILS_H

#include "tsp.h"

double second();
double rand01();

int compar(const void *a, const void *b);
int is_exact_method(const char *method);
int runs(int argc, char **argv);

void copy_sol(const solution *src, solution *dst);
void set_seed(int seed);
void update_perf_csv(const instance *inst, double *run_costs, int num_runs);

#endif /* UTILS_H */
