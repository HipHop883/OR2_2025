#include "tsp.h"
#include "tsp_cplex.h"
#include "tsp_greedy.h"
#include "utils.h"
#include <cplex.h>

/**
 * Calculates the index of the variable x(i, j) in the CPLEX model
 * @param i index of node i
 * @param j index of node j
 * @param inst instance
 * @return index of the variable x(i, j) to be used in the xstar array
 */
int xpos(int i, int j, const instance *inst)
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

int add_warm_start(CPXENVptr env, CPXLPptr lp, const instance *inst, const solution *sol)
{
    int n = inst->nnodes;
    int ncols = CPXgetnumcols(env, lp);
    int *indices = malloc(ncols * sizeof(int));
    double *values = calloc(ncols, sizeof(double));

    if (!indices || !values)
    {
        print_error("Memory allocation failed in add_warm_start");
        return EXIT_FAILURE;
    }

    for (int i = 0; i < n; i++)
    {
        int j = sol->tour[i + 1];
        if (i != j)
            values[xpos(i, j, inst)] = 1.0;
    }

    for (int i = 0; i < ncols; i++)
        indices[i] = i;

    int beg = 0;
    int effort = CPX_MIPSTART_AUTO;
    char *name = "heuristic_start";

    int status = CPXaddmipstarts(env, lp, 1, ncols, &beg, indices, values, &effort, &name);
    if (status)
    {
        print_error("CPXaddmipstarts failed");
        return EXIT_FAILURE;
    }

    free(indices);
    free(values);
    return EXIT_SUCCESS;
}

int patch_solution(const double *xstar, instance *inst, solution *sol, int *comp, int ncomp)
{
    int n = inst->nnodes;

    int **adj = malloc(n * sizeof(int *));
    int *degree = calloc(n, sizeof(int));
    for (int i = 0; i < n; i++)
    {
        adj[i] = malloc(2 * sizeof(int));
        degree[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (xstar[xpos(i, j, inst)] > 0.5 + EPSILON)
            {
                adj[i][degree[i]++] = j;
                adj[j][degree[j]++] = i;
            }
        }
    }

    double best_extra_cost = CPX_INFBOUND;
    int best_i = -1, best_j = -1, best_u = -1, best_v = -1;

    for (int i = 0; i < n; i++)
    {
        if (degree[i] != 2)
            continue;
        int j = adj[i][0];
        if (comp[i] != comp[j])
            continue;

        for (int u = 0; u < n; u++)
        {
            if (u == i || u == j || comp[u] == comp[i])
                continue;

            for (int v = 0; v < n; v++)
            {
                if (v == i || v == j || v == u || comp[v] == comp[i] || comp[v] == comp[u])
                    continue;

                double removed_cost = inst->cost_matrix[i][j];
                double added_cost = inst->cost_matrix[i][u] + inst->cost_matrix[j][v];
                double extra = added_cost - removed_cost;

                if (extra < best_extra_cost)
                {
                    best_extra_cost = extra;
                    best_i = i;
                    best_j = j;
                    best_u = u;
                    best_v = v;
                }
            }
        }
    }

    if (best_i == -1)
    {
        printf("No patch found — fallback heuristic failed\n");
        for (int i = 0; i < n; i++)
            free(adj[i]);
        free(adj);
        free(degree);
        return EXIT_FAILURE;
    }

    printf("Patching: removing (%d,%d), adding (%d,%d) and (%d,%d) with extra cost %.2f\n",
           best_i, best_j, best_i, best_u, best_j, best_v, best_extra_cost);

    for (int d = 0; d < 2; d++)
    {
        if (adj[best_i][d] == best_j)
        {
            adj[best_i][d] = best_u;
            break;
        }
    }
    for (int d = 0; d < 2; d++)
    {
        if (adj[best_j][d] == best_i)
        {
            adj[best_j][d] = best_v;
            break;
        }
    }

    adj[best_u][degree[best_u]++] = best_i;
    adj[best_v][degree[best_v]++] = best_j;

    int *tour = malloc((n + 1) * sizeof(int));
    int *visited = calloc(n, sizeof(int));
    int current = 0;
    tour[0] = current;
    visited[current] = 1;

    for (int k = 1; k < n; k++)
    {
        int next = (visited[adj[current][0]] == 0) ? adj[current][0] : adj[current][1];
        tour[k] = next;
        visited[next] = 1;
        current = next;
    }
    tour[n] = tour[0];

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
 * Callback function for CPLEX to add SEC constraints
 * @param context CPLEX callback context
 * @param contextid context ID
 * @param userhandle user data
 *
 * @return 0 if successful, 1 otherwise
 */
static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle)
{
    instance *inst = (instance *)userhandle;

    int n = inst->nnodes;
    int max_edges = n * (n - 1) / 2;

    double *xstar = (double *)malloc(inst->ncols * sizeof(double));
    if (!xstar)
    {
        print_error("malloc failed for xstar");
        return 1;
    }

    double objval = CPX_INFBOUND;
    if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval))
    {
        print_error("CPXcallbackgetcandidatepoint error");
        free(xstar);
        return 1;
    }

    int *comp = (int *)malloc(n * sizeof(int));
    if (!comp)
    {
        print_error("malloc failed for comp");
        free(xstar);
        return 1;
    }

    int ncomp = build_components(xstar, comp, inst);
    if (ncomp <= 1)
    {
        free(comp);
        free(xstar);
        return 0; // feasible
    }

    // Add one SEC per component
    int izero = 0;
    char sense = 'L';

    for (int k = 0; k < ncomp; k++)
    {
        int *index = (int *)malloc(max_edges * sizeof(int));
        double *value = (double *)malloc(max_edges * sizeof(double));
        if (!index || !value)
        {
            print_error("malloc failed for cut data");
            free(comp);
            free(xstar);
            return 1;
        }

        int nnz = 0;
        int rhs = -1;

        for (int i = 0; i < n; i++)
        {
            if (comp[i] != k)
                continue;
            rhs++;

            for (int j = i + 1; j < n; j++)
            {
                if (comp[j] != k)
                    continue;
                index[nnz] = xpos(i, j, inst);
                value[nnz++] = 1.0;

                // Safety check
                if (nnz >= max_edges)
                {
                    print_error("SEC constraint too large — aborting");
                    free(index);
                    free(value);
                    free(comp);
                    free(xstar);
                    return 1;
                }
            }
        }

        if (rhs >= 1 && nnz > 0)
        {
            double rhs_val = (double)rhs;
            if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs_val, &sense, &izero, index, value))
            {
                print_error("CPXcallbackrejectcandidate error");
                free(index);
                free(value);
                free(comp);
                free(xstar);
                return 1;
            }
        }

        free(index);
        free(value);
    }

    free(comp);
    free(xstar);
    return 0;
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

    inst->ncols = CPXgetnumcols(env, lp);

    solution warm_sol;
    warm_sol.tour = (int *)malloc((inst->nnodes + 1) * sizeof(int));
    if (warm_sol.tour == NULL)
    {
        print_error("Memory allocation failed for warm_sol->tour");
        return EXIT_FAILURE;
    }
    warm_sol.initialized = 0;

    if (apply_greedy_search(inst, &warm_sol) == 0)
    {
        if (add_warm_start(env, lp, inst, &warm_sol))
            printf("Warning: failed to add warm start\n");
        else if (VERBOSE >= 50)
            printf("Warm start injected with cost %.2f\n", warm_sol.cost);
        free_sol(&warm_sol);
    }
    else
    {
        printf("Warning: greedy warm start failed\n");
    }

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

            if (ncomp > 1)
            {
                if (patch_solution(xstar, inst, sol, comp, ncomp))
                    return EXIT_FAILURE;
            }

            break;
        }

        CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time);
        CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
        CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);
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

/**
 * Apply CPLEX to solve the TSP instance using Branch and Cut method
 * @param inst instance
 * @param sol solution
 *
 * @return 0 if successful, 1 otherwise
 */
int apply_cplex_branchcut(instance *inst, solution *sol)
{
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
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    if (build_model(inst, env, lp))
        return EXIT_FAILURE;

    // Save ncols to use in callback
    inst->ncols = CPXgetnumcols(env, lp);

    // Install lazy constraint callback (CANDIDATE context)
    CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
    if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst))
    {
        print_error("CPXcallbacksetfunc error");
        return EXIT_FAILURE;
    }

    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);
    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
    CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);

    if (CPXmipopt(env, lp))
    {
        print_error("CPXmipopt error");
        return EXIT_FAILURE;
    }

    // Extract final solution
    double *xstar = calloc(inst->ncols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, inst->ncols - 1))
    {
        print_error("CPXgetx error");
        return EXIT_FAILURE;
    }

    if (build_solution(xstar, inst, sol))
    {
        print_error("Failed to build solution");
        return EXIT_FAILURE;
    }

    free(xstar);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return EXIT_SUCCESS;
}
