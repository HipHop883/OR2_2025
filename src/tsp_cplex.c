#include <cplex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tsp.h"
#include "tsp_cplex.h"

#define EPS 1e-5

int xpos(int i, int j, instance *inst)
{
    if (i == j)
        print_error("xpos called with i == j");
    if (i > j)
        return xpos(j, i, inst);
    return i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
}

void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
{
    int n = inst->nnodes;
    int izero = 0;
    char binary = 'B';

    char **cname = (char **)calloc(1, sizeof(char *));
    cname[0] = (char *)calloc(100, sizeof(char));

    // Add binary vars x(i,j) for i < j
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
            double obj = inst->cost_matrix[i][j];
            double lb = 0.0;
            double ub = 1.0;

            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
                print_error("CPXnewcols failed on x var.s");

            if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst))
                print_error("Wrong index in xpos()");
        }
    }

    // Degree constraints
    int *index = (int *)malloc(n * sizeof(int));
    double *value = (double *)malloc(n * sizeof(double));

    for (int h = 0; h < n; h++)
    {
        double rhs = 2.0;
        char sense = 'E';
        sprintf(cname[0], "degree(%d)", h + 1);
        int nnz = 0;

        for (int i = 0; i < n; i++)
        {
            if (i == h)
                continue;
            index[nnz] = xpos(i, h, inst);
            value[nnz++] = 1.0;
        }

        if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname))
            print_error("CPXaddrows failed [degree]");
    }

    free(value);
    free(index);
    free(cname[0]);
    free(cname);

    if (VERBOSE >= 60)
        CPXwriteprob(env, lp, "model.lp", NULL);
}

void build_solution(const double *xstar, instance *inst, solution *sol)
{
    int n = inst->nnodes;
    int *next = (int *)malloc(n * sizeof(int));
    int *visited = (int *)calloc(n, sizeof(int));

    for (int i = 0; i < n; i++)
        next[i] = -1;

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            int k = xpos(i, j, inst);
            if (xstar[k] > 0.5)
            {
                if (next[i] == -1)
                    next[i] = j;
                else
                    next[j] = i;
            }
        }
    }

    int *tour = (int *)malloc((n + 1) * sizeof(int));
    int curr = 0, count = 0;
    while (count < n)
    {
        tour[count++] = curr;
        visited[curr] = 1;
        curr = next[curr];
        if (visited[curr])
            break;
    }
    tour[n] = tour[0]; // close the tour

    sol->tour = tour;
    sol->cost = 0.0;
    for (int i = 0; i < n; i++)
        sol->cost += inst->cost_matrix[tour[i]][tour[i + 1]];

    sol->initialized = 1;

    free(next);
    free(visited);
}

int apply_cplex(instance *inst, solution *sol)
{
    int error;

    CPXENVptr env = CPXopenCPLEX(&error);
    if (env == NULL)
        print_error("CPXopenCPLEX error");

    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
    if (lp == NULL)
        print_error("CPXcreateprob error");

    build_model(inst, env, lp);

    // Parameters
    CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);
    CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);
    CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 2);
    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);

    // Solve
    if (CPXmipopt(env, lp))
        print_error("CPXmipopt error");

    int ncols = CPXgetnumcols(env, lp);
    double *xstar = (double *)calloc(ncols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, ncols - 1))
        print_error("CPXgetx error");

    build_solution(xstar, inst, sol);

    free(xstar);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return 0;
}
