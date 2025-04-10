#include "tsp_cplex.h"
#include "utils.h"
#include <cplex.h>

#define EPS 1e-5

/**
 * Calculates the index of the variable x(i, j) in the CPLEX model
 * @param i index of node i
 * @param j index of node j
 * @param inst instance
 * @return index of the variable x(i, j) to be used in the xstar array
 */
int xpos(int i, int j, instance *inst)
{
    if (i == j)
        print_error("xpos called with i == j");
    if (i > j)
        return xpos(j, i, inst);
    return i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
}

/**
 * Builds the CPLEX model
 * @param inst instance
 * @param env CPLEX environment
 * @param lp CPLEX problem
 * @return 0 if successful, 1 otherwise
 */
int build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
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
            {
                print_error("CPXnewcols failed on x var.s");
                return EXIT_FAILURE;
            }

            if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst))
            {
                print_error("Wrong index in xpos()");
                return EXIT_FAILURE;
            }
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
        {
            print_error("CPXaddrows failed [degree]");
            return EXIT_FAILURE;
        }
    }

    free(value);
    free(index);
    free(cname[0]);
    free(cname);

    if (VERBOSE >= 60)
        CPXwriteprob(env, lp, "model.lp", NULL);

    return EXIT_SUCCESS;
}

/**
 * Builds the solution from the xstar vector
 * @param xstar solution vector
 * @param inst instance
 * @param sol solution to be filled
 * @return 0 if successful, 1 otherwise
 */
int build_solution(const double *xstar, instance *inst, solution *sol)
{
    int n = inst->nnodes;

    // Build adjacency list: each node must have exactly 2 neighbors
    int **adj = malloc(n * sizeof(int *));
    int *degree = calloc(n, sizeof(int));
    for (int i = 0; i < n; i++)
    {
        adj[i] = malloc(2 * sizeof(int)); // Each node must have 2 neighbors
        degree[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (xstar[xpos(i, j, inst)] > 0.5 + EPSILON)
            {
                if (degree[i] >= 2 || degree[j] >= 2)
                {
                    print_error("Node degree exceeds 2 — invalid solution");
                    return EXIT_FAILURE;
                }
                adj[i][degree[i]++] = j;
                adj[j][degree[j]++] = i;
            }
        }
    }

    // Validate all degrees are exactly 2
    for (int i = 0; i < n; i++)
    {
        if (degree[i] != 2)
        {
            printf("Node %d has degree %d instead of 2 - not a valid tour\n", i, degree[i]);
            return EXIT_FAILURE;
        }
    }

    // Reconstruct tour
    int *tour = malloc((n + 1) * sizeof(int));
    int *visited = calloc(n, sizeof(int));
    int curr = 0;
    tour[0] = curr;
    visited[curr] = 1;

    // 
    for (int i = 1; i < n; i++)
    {
        int next = (visited[adj[curr][0]] == 0) ? adj[curr][0] : adj[curr][1];
        tour[i] = next;
        visited[next] = 1;
        curr = next;
    }

    tour[n] = tour[0]; // close tour

    free_sol(sol);
    sol->tour = tour;
    sol->initialized = 1;
    evaluate_path_cost(inst, sol);

    for (int i = 0; i < n; i++)
        free(adj[i]);
    free(adj);
    free(degree);
    free(visited);

    return EXIT_SUCCESS;
}

/**
 * Add SEC constraints to the model
 * @param env CPLEX environment
 * @param lp CPLEX problem
 * @param comp component array
 * @param ncomp number of components
 * @param inst instance
 * @return 0 if successful, 1 otherwise 
 */
int add_sec(CPXENVptr env, CPXLPptr lp, int *comp, int ncomp, instance *inst)
{
    int n = inst->nnodes;
    int izero = 0;
    char sense = 'L';

    for (int k = 0; k < ncomp; k++)
    {
        // Upper bound estimate on SEC size
        int max_nz = n * n;
        int *index = (int *)malloc(max_nz * sizeof(int));
        double *value = (double *)malloc(max_nz * sizeof(double));
        if (!index || !value)
        {
            print_error("Memory allocation failed in add_sec");
            return EXIT_FAILURE;
        }

        int nnz = 0;
        int rhs = -1;

        for (int i = 0; i < n; i++)
        {
            if (comp[i] != k)
                continue;
            rhs++;

            for (int j = 0; j < n; j++)
            {
                if (i == j)
                    continue;
                if (comp[j] != k)
                    continue;

                index[nnz] = xpos(i, j, inst);
                value[nnz++] = 1.0;
            }
        }

        if (rhs >= 1 && nnz > 0)
        {
            int status = CPXaddrows(env, lp,
                                    0,   // no new rows in A matrix
                                    1,   // 1 new constraint
                                    nnz, // number of non-zero elements
                                    (double[]){(double)rhs},
                                    &sense,
                                    &izero,
                                    index,
                                    value,
                                    NULL,
                                    NULL);

            if (status)
            {
                print_error("CPXaddrows failed in add_sec");
                return EXIT_FAILURE;
            }

            if (VERBOSE >= 50)
            {
                printf(" SEC added for component %d with %d nodes, %d edges, RHS = %d\n",
                       k, rhs + 1, nnz, rhs);
            }
        }

        free(index);
        free(value);
    }

    return EXIT_SUCCESS;
}

/**
 * Assigns components to each node based on the solution vector xstar using the BFS algorithm
 * @param xstar solution vector
 * @param comp component array to be filled
 * @param inst instance
 * @return number of components found
 */
int build_components(const double *xstar, int *comp, instance *inst)
{
    int n = inst->nnodes;
    int ncomp = 0;

    // Initialize all nodes as unvisited
    for (int i = 0; i < n; i++)
        comp[i] = -1;

    int *queue = (int *)malloc(n * sizeof(int));
    if (!queue)
    {
        print_error("Memory allocation failed in build_components");
        return EXIT_FAILURE;
    }

    for (int start = 0; start < n; start++)
    {
        // Skip if already visited
        if (comp[start] != -1)
            continue;

        int q_head = 0, q_tail = 0;
        queue[q_tail++] = start;
        comp[start] = ncomp;

        while (q_head < q_tail)
        {
            int curr = queue[q_head++];

            for (int j = 0; j < n; j++)
            {
                if (j == curr)
                    continue;

                if (comp[j] == -1 && xstar[xpos(curr, j, inst)] > 0.5)
                {
                    comp[j] = ncomp;
                    queue[q_tail++] = j;
                }
            }
        }

        ncomp++;
    }

    free(queue);
    return ncomp;
}

/**
 * Apply CPLEX to solve the TSP instance using the Benders decomposition method
 * @param inst instance
 * @param sol solution
 * @return 0 if successful, 1 otherwise
 */
int apply_cplex_beneders(instance *inst, solution *sol)
{
    int n = inst->nnodes;
    double start_time = second();

    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    if (env == NULL)
    {
        print_error("CPXopenCPLEX error");
        return EXIT_FAILURE;
    }

    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
    if (lp == NULL)
    {
        print_error("CPXcreateprob error");
        return EXIT_FAILURE;
    }

    if (build_model(inst, env, lp))
        return EXIT_FAILURE;

    int *comp = malloc(n * sizeof(int));
    double *xstar = NULL;
    int ncomp = -1;
    int iteration = 0;

    while (1)
    {
        iteration++;

        double elapsed = second() - start_time;
        double remaining_time = inst->timelimit - elapsed;
        if (remaining_time <= 0.0)
        {
            print_error("Time limit reached before full tour found.");
            break;
        }

        CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time);
        CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 2); // optional
        if (CPXmipopt(env, lp))
        {
            print_error("CPXmipopt error");
            return EXIT_FAILURE;
        }

        int ncols = CPXgetnumcols(env, lp);
        xstar = (double *)calloc(ncols, sizeof(double));
        if (CPXgetx(env, lp, xstar, 0, ncols - 1))
        {
            print_error("CPXgetx error");
            return EXIT_FAILURE;
        }

        ncomp = build_components(xstar, comp, inst);

        printf("Iteration %d | Components = %d | Time = %.2lf sec\n", iteration, ncomp, second() - start_time);

        if (ncomp == 1)
        {
            if (build_solution(xstar, inst, sol))
                return EXIT_FAILURE;
            break;
        }

        if (add_sec(env, lp, comp, ncomp, inst))
            return EXIT_FAILURE;
        free(xstar);
    }

    free(comp);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return EXIT_SUCCESS;
}
