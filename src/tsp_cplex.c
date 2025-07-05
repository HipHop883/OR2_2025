#include "tsp.h"
#include "tsp_cplex.h"
#include "tsp_greedy.h"
#include "utils.h"
#include <assert.h>
#include <cplex.h>

#define ASSERT_INDEX_IN_RANGE(idx, max) assert((idx) >= 0 && (idx) < (max))
#define MAX_DEGREE 2

#ifndef MAX_MIN_MACROS
#define MAX_MIN_MACROS
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

/**
 * Calculates the index of the variable x(i, j) in the CPLEX model
 * @param i index of node i
 * @param j index of node j
 * @param inst instance
 * @return index of the variable x(i, j) to be used in the xstar array
 */
static inline int xpos(int i, int j, const instance *inst)
{
    assert(inst != NULL);
    assert(i != j);

    ASSERT_INDEX_IN_RANGE(i, inst->nnodes);
    ASSERT_INDEX_IN_RANGE(j, inst->nnodes);

    if (i > j)
    {
        int tmp = i;
        i = j;
        j = tmp;
    }

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
    const int n = inst->nnodes;
    int izero = 0;
    char binary = 'B';
    char cname[100];
    char *nameptr = cname; // pointer for passing to CPLEX API

    // === Add binary variables for x(i,j), i < j ===
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            snprintf(cname, sizeof(cname), "x(%d,%d)", i + 1, j + 1);
            double obj = inst->cost_matrix[i][j];
            double lb = 0.0;
            double ub = 1.0;

            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, &nameptr))
            {
                print_error("CPXnewcols failed");
                return EXIT_FAILURE;
            }

            if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst))
            {
                print_error("Mismatch in xpos index");
                return EXIT_FAILURE;
            }
        }
    }

    // === Add degree constraints ===
    int *index = malloc(n * sizeof(int));
    double *value = malloc(n * sizeof(double));

    if (!index || !value)
    {
        print_error("Memory allocation failed");
        return EXIT_FAILURE;
    }

    for (int h = 0; h < n; h++)
    {
        double rhs = 2.0;
        char sense = 'E';
        snprintf(cname, sizeof(cname), "degree(%d)", h + 1);
        int nnz = 0;

        for (int i = 0; i < n; i++)
        {
            if (i == h)
                continue;
            index[nnz] = xpos(i, h, inst);
            value[nnz++] = 1.0;
        }

        if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &nameptr))
        {
            print_error("CPXaddrows failed [degree constraint]");
            free(index);
            free(value);
            return EXIT_FAILURE;
        }
    }

    if (VERBOSE >= 50)
        printf("[CPLEX] Model built with %d variables for %d nodes\n", CPXgetnumcols(env, lp), n);

    free(index);
    free(value);

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
    int **adj = malloc(n * sizeof(int *));
    int *degree = calloc(n, sizeof(int));

    if (!adj || !degree)
    {
        print_error("Memory allocation failed [adj or degree]");
        return EXIT_FAILURE;
    }

    for (int i = 0; i < n; i++)
    {
        adj[i] = malloc(2 * sizeof(int));
        if (!adj[i])
        {
            print_error("Memory allocation failed [adj[i]]");
            return EXIT_FAILURE;
        }
        degree[i] = 0;
    }

    // === Build adjacency list from xstar ===
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (xstar[xpos(i, j, inst)] > 0.5 + EPSILON)
            {
                if (degree[i] >= 2 || degree[j] >= 2)
                {
                    print_error("Invalid degree — exceeds 2");
                    return EXIT_FAILURE;
                }
                adj[i][degree[i]++] = j;
                adj[j][degree[j]++] = i;
            }
        }
    }

    // === Validate all degrees ===
    for (int i = 0; i < n; i++)
    {
        if (degree[i] != 2)
        {
            printf("Node %d has degree %d instead of 2 - not a valid tour\n", i, degree[i]);
            return EXIT_FAILURE;
        }
    }

    // === Reconstruct tour ===
    int *tour = malloc((n + 1) * sizeof(int));
    int *visited = calloc(n, sizeof(int));
    if (!tour || !visited)
    {
        print_error("Memory allocation failed [tour or visited]");
        return EXIT_FAILURE;
    }

    int curr = 0;
    tour[0] = curr;
    visited[curr] = 1;

    for (int i = 1; i < n; i++)
    {
        int next = (visited[adj[curr][0]] == 0) ? adj[curr][0] : adj[curr][1];
        tour[i] = next;
        visited[next] = 1;
        curr = next;
    }
    tour[n] = tour[0];

    if (VERBOSE >= 50)
        printf("[CPLEX] Tour reconstructed successfully from xstar\n");

    // === Finalize ===
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
        int *index = malloc(n * n * sizeof(int));
        double *value = malloc(n * n * sizeof(double));
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

            for (int j = i + 1; j < n; j++)
            {
                if (comp[j] != k)
                    continue;

                index[nnz] = xpos(i, j, inst);
                value[nnz++] = 1.0;

                if (nnz >= inst->ncols)
                {
                    print_error("SEC too large — too many edges");
                    free(index);
                    free(value);
                    return EXIT_FAILURE;
                }
            }
        }

        if (rhs >= 1 && nnz > 0)
        {
            double rhs_val = (double)rhs;
            if (CPXaddrows(env, lp, 0, 1, nnz, &rhs_val, &sense, &izero, index, value, NULL, NULL))
            {
                print_error("CPXaddrows failed in add_sec");
                free(index);
                free(value);
                return EXIT_FAILURE;
            }

            if (VERBOSE >= 50)
                printf("[CPLEX] SEC added for comp %d: %d nodes, %d terms, rhs = %d\n", k, rhs + 1, nnz, rhs);
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

    int *queue = malloc(n * sizeof(int));
    if (!queue)
    {
        print_error("Memory allocation failed in build_components");
        return EXIT_FAILURE;
    }

    for (int start = 0; start < n; start++)
    {
        if (comp[start] != -1)
            continue;

        // Start a new BFS for a new component
        int q_head = 0, q_tail = 0;
        queue[q_tail++] = start;
        comp[start] = ncomp;

        while (q_head < q_tail)
        {
            int curr = queue[q_head++];
            for (int j = 0; j < n; j++)
            {
                if (j == curr || comp[j] != -1)
                    continue;
                if (xstar[xpos(curr, j, inst)] > 0.5 + EPSILON)
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
 * Adds a warm start to the CPLEX model using the provided solution
 *
 * @param env CPLEX environment
 * @param lp CPLEX problem
 * @param inst instance
 * @param sol solution to be used for warm start
 * @return 0 if successful, 1 otherwise
 */
int add_warm_start(CPXENVptr env, CPXLPptr lp, const instance *inst, const solution *sol)
{
    if (!sol || !sol->initialized || !sol->tour)
    {
        print_error("Invalid or uninitialized solution in add_warm_start");
        return EXIT_FAILURE;
    }

    int n = inst->nnodes;
    int ncols = CPXgetnumcols(env, lp);

    if (ncols <= 0)
    {
        print_error("CPXgetnumcols returned invalid column count");
        return EXIT_FAILURE;
    }

    int *indices = malloc(ncols * sizeof(int));
    double *values = calloc(ncols, sizeof(double));

    if (!indices || !values)
    {
        print_error("Memory allocation failed in add_warm_start");
        return EXIT_FAILURE;
    }

    // Set all indices
    for (int i = 0; i < ncols; ++i)
        indices[i] = i;

    // Convert the tour into x variable values
    for (int i = 0; i < n; i++)
    {
        int from = sol->tour[i];
        int to = sol->tour[i + 1];

        if (from < 0 || from >= n || to < 0 || to >= n || from == to)
        {
            fprintf(stderr, "Invalid tour edge: %d → %d (n = %d)\n", from, to, n);
            continue;
        }

        values[xpos(from, to, inst)] = 1.0;
    }

    int beg = 0;
    int effort = CPX_MIPSTART_AUTO;
    char *mipstart_name = "heuristic_start";

    int status = CPXaddmipstarts(env, lp, 1, ncols, &beg, indices, values, &effort, &mipstart_name);
    if (status)
    {
        fprintf(stderr, "Warning: CPXaddmipstarts failed (status %d)\n", status);
        free(indices);
        free(values);
        return EXIT_FAILURE;
    }
    else if (VERBOSE >= 50)
        printf("[CPLEX] Warm start added successfully with cost %.2f\n", sol->cost);

    free(indices);
    free(values);
    return EXIT_SUCCESS;
}

/**
 * Attempts to patch a non-Hamiltonian solution into a feasible tour by
 * reconnecting disconnected components.
 *
 * @param xstar solution vector
 * @param inst problem instance
 * @param sol solution to be constructed
 * @param comp array indicating component index of each node
 * @param ncomp number of components
 *
 * @return EXIT_SUCCESS if successful, EXIT_FAILURE otherwise
 */
int patch_solution(const double *xstar_in, instance *inst, solution *sol)
{
    int n = inst->nnodes;

    // === Allocate support data structures ===
    double *xstar = malloc(inst->ncols * sizeof(double));
    int *comp = malloc(n * sizeof(int));

    int **adj = malloc(n * sizeof(int *));
    int *degree = calloc(n, sizeof(int));
    if (!adj || !degree)
    {
        print_error("Memory allocation failed [adj or degree]");
        free(adj);
        free(degree);
        return EXIT_FAILURE;
    }

    int *visited;
    int *tour;
    int did_patch = 0;

    memcpy(xstar, xstar_in, inst->ncols * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        adj[i] = malloc(MAX_DEGREE * sizeof(int));
        if (!adj[i])
        {
            print_error("Memory allocation failed [adj[i]]");
            for (int j = 0; j < i; j++)
                free(adj[j]);
            free(adj);
            free(degree);
            return EXIT_FAILURE;
        }
        degree[i] = 0;
    }

    // === Build adjacency from xstar ===
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

    int ncomp = build_components(xstar, comp, inst);

    while (ncomp > 1)
    {
        double best_extra_cost = CPX_INFBOUND;
        int bi = -1, bj = -1, bu = -1, bv = -1;

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (comp[i] == comp[j])
                    continue;

                for (int di = 0; di < degree[i]; di++)
                {
                    int u = adj[i][di];

                    for (int dj = 0; dj < degree[j]; dj++)
                    {
                        int v = adj[j][dj];
                        if (u == j && v == i)
                            continue;
                        if (u == v)
                            continue;
                        double extra_cost = inst->cost_matrix[i][j] + inst->cost_matrix[u][v] - inst->cost_matrix[i][u] - inst->cost_matrix[j][v];
                        if (extra_cost < best_extra_cost)
                        {
                            best_extra_cost = extra_cost;
                            bi = i;
                            bj = j;
                            bu = u;
                            bv = v;
                        }
                    }
                }
            }
        }

        if (bi < 0)
        {
            if (VERBOSE >= 20)
                printf("[PATCH] No patch found — fallback heuristic failed\n");

            for (int i = 0; i < n; i++)
                free(adj[i]);
            free(adj);
            free(degree);
            return EXIT_FAILURE;
        }

        if (VERBOSE >= 30)
            printf("[PATCH] Removing (%d,%d), adding (%d,%d) and (%d,%d) | Delta cost = %.2f\n", bi, bj, bi, bu, bj, bv, best_extra_cost);

        // === Apply patch ===
        for (int k = 0; k < degree[bi]; k++)
        {
            if (adj[bi][k] == bu)
            {
                adj[bi][k] = adj[bi][--degree[bi]];
                break;
            }
        }
        for (int k = 0; k < degree[bu]; k++)
        {
            if (adj[bu][k] == bi)
            {
                adj[bu][k] = adj[bu][--degree[bu]];
                break;
            }
        }

        for (int k = 0; k < degree[bj]; k++)
        {
            if (adj[bj][k] == bv)
            {
                adj[bj][k] = adj[bj][--degree[bj]];
                break;
            }
        }
        for (int k = 0; k < degree[bv]; k++)
        {
            if (adj[bv][k] == bj)
            {
                adj[bv][k] = adj[bv][--degree[bv]];
                break;
            }
        }

        adj[bi][degree[bi]++] = bj;
        adj[bj][degree[bj]++] = bi;
        adj[bu][degree[bu]++] = bv;
        adj[bv][degree[bv]++] = bu;

        did_patch = 1;

        // === Rebuild xstar from adj ===
        memset(xstar, 0, inst->ncols * sizeof(double));
        for (int i = 0; i < n; i++)
        {
            for (int d = 0; d < degree[i]; d++)
            {
                int j = adj[i][d];
                if (i < j)
                    xstar[xpos(i, j, inst)] = 1.0;
            }
        }

        // === Compute new ncomp ===
        ncomp = build_components(xstar, comp, inst);
    }

    if (!did_patch)
    {
        print_error("No valid patch found");
        for (int i = 0; i < n; i++)
            free(adj[i]);
        free(adj);
        free(degree);
        free(comp);
        free(xstar);

        return EXIT_FAILURE;
    }

    // === Reconstruct tour ===
    tour = malloc((n + 1) * sizeof(int));
    visited = calloc(n, sizeof(int));
    int curr = 0;
    tour[0] = curr;
    visited[curr] = 1;

    for (int k = 1; k < n; k++)
    {
        int a = adj[curr][0], b = adj[curr][1];
        int next = visited[a] ? b : a;
        tour[k] = next;
        visited[next] = 1;
        curr = next;
    }
    tour[n] = tour[0];

    // === Finalize solution ===
    free_sol(sol);
    sol->tour = tour;
    sol->initialized = 1;
    evaluate_path_cost(inst, sol);
    for (int i = 0; i < n; i++)
        free(adj[i]);
    free(adj);
    free(degree);
    free(comp);
    free(xstar);
    if (visited)
        free(visited);
    return (sol->initialized ? EXIT_SUCCESS : EXIT_FAILURE);
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
    const int n = inst->nnodes;
    const int max_edges = n * (n - 1) / 2;
    int izero = 0;
    char sense = 'L';

    // === Allocate memory for solution vector ===
    double *xstar = malloc(inst->ncols * sizeof(double));
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

    int *comp = malloc(n * sizeof(int));
    if (!comp)
    {
        print_error("malloc failed for comp");
        free(xstar);
        return 1;
    }

    int ncomp = build_components(xstar, comp, inst);
    if (ncomp <= 1)
    {
        // Feasible — no subtours
        free(comp);
        free(xstar);
        return 0;
    }

    // === Loop over components to add lazy SEC constraints ===
    for (int k = 0; k < ncomp; k++)
    {
        int *index = malloc(max_edges * sizeof(int));
        double *value = malloc(max_edges * sizeof(double));
        if (!index || !value)
        {
            print_error("malloc failed for SEC cut arrays");
            free(index);
            free(value);
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

                if (nnz >= max_edges)
                {
                    print_error("SEC too large — aborting");
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
int apply_cplex_benders(instance *inst, solution *sol)
{
    const int n = inst->nnodes;
    const double start_time = second();

    // === Initialize CPLEX ===
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    if (!env)
    {
        print_error("CPXopenCPLEX error");
        return EXIT_FAILURE;
    }

    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
    if (!lp)
    {
        print_error("CPXcreateprob error");
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // === Build initial model ===
    if (build_model(inst, env, lp))
    {
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    inst->ncols = CPXgetnumcols(env, lp);

    // === Warm start ===
    if (sol->initialized && sol->tour)
    {
        if (add_warm_start(env, lp, inst, sol))
            fprintf(stderr, "Warning: failed to add warm start\n");
        else if (VERBOSE >= 50)
            printf("Warm start injected with cost %.2f\n", sol->cost);
    }
    else
    {
        if (VERBOSE >= 50)
            printf("No valid warm start solution provided.\n");
    }

    int *comp = malloc(n * sizeof(int));
    if (!comp)
    {
        print_error("Memory allocation failed [comp]");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    double *xstar = NULL;
    int ncomp = -1;
    int iteration = 0;

    // === Main optimization loop ===
    while (1)
    {
        iteration++;

        double elapsed = second() - start_time;
        double remaining_time = inst->timelimit - elapsed;
        if (remaining_time <= 0.0)
        {
            fprintf(stderr, "Time limit reached before full tour found.\n");

            if (ncomp > 1 && xstar)
            {
                if (patch_solution(xstar, inst, sol))
                {
                    free(comp);
                    free(xstar);
                    CPXfreeprob(env, &lp);
                    CPXcloseCPLEX(&env);
                    return EXIT_FAILURE;
                }
            }

            break;
        }

        CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time);
        CPXsetintparam(env, CPX_PARAM_SCRIND, VERBOSE >= 100 ? CPX_ON : CPX_OFF); // Turn on\off CPLEX output for iterations
        CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, VERBOSE >= 100 ? 4 : 0);        // Minimal display

        if (CPXmipopt(env, lp))
        {
            print_error("CPXmipopt error");
            free(comp);
            free(xstar);
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return EXIT_FAILURE;
        }

        int ncols = CPXgetnumcols(env, lp);

        if (xstar)
            free(xstar);
        xstar = calloc(ncols, sizeof(double));

        if (!xstar)
        {
            print_error("Memory allocation failed [xstar]");
            free(comp);
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return EXIT_FAILURE;
        }

        if (CPXgetx(env, lp, xstar, 0, ncols - 1))
        {
            print_error("CPXgetx error");
            free(comp);
            free(xstar);
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return EXIT_FAILURE;
        }

        ncomp = build_components(xstar, comp, inst);

        if (VERBOSE >= 50)
            printf("[BENDERS] Iteration %2d | Components: %2d | Elapsed: %.2fs\n", iteration, ncomp, second() - start_time);

        if (ncomp == 1)
        {
            if (build_solution(xstar, inst, sol))
            {
                free(comp);
                free(xstar);
                CPXfreeprob(env, &lp);
                CPXcloseCPLEX(&env);
                return EXIT_FAILURE;
            }
            break;
        }

        if (add_sec(env, lp, comp, ncomp, inst))
        {
            free(comp);
            free(xstar);
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return EXIT_FAILURE;
        }
    }

    if (VERBOSE >= 10)
        printf("[BENDERS] Optimization finished. Final cost: %.2f\n", sol->cost);

    // === Cleanup ===
    free(comp);
    free(xstar);
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
    // === Initialize CPLEX ===
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    if (!env)
    {
        print_error("CPXopenCPLEX error");
        return EXIT_FAILURE;
    }

    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
    if (!lp)
    {
        print_error("CPXcreateprob error");
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // === Build model ===
    if (build_model(inst, env, lp))
    {
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    inst->ncols = CPXgetnumcols(env, lp);

    // === Set up callback for lazy constraints (SEC) ===
    CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
    if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst))
    {
        print_error("CPXcallbacksetfunc error");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // === Configure CPLEX parameters ===
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);
    CPXsetintparam(env, CPX_PARAM_SCRIND, VERBOSE >= 100 ? CPX_ON : CPX_OFF); // Turn on\off CPLEX output for iterations
    CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, VERBOSE >= 100 ? 4 : 0);        // Minimal display

    // === Solve the model ===
    if (CPXmipopt(env, lp))
    {
        print_error("CPXmipopt error");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // === Extract solution ===
    double *xstar = calloc(inst->ncols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, inst->ncols - 1))
    {
        print_error("CPXgetx error");
        free(xstar);
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    if (build_solution(xstar, inst, sol))
    {
        print_error("Failed to build solution");
        free(xstar);
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    if (VERBOSE >= 10)
        printf("[BRANCHCUT] Optimization complete. Final cost: %.2f\n", sol->cost);

    // === Cleanup ===
    free(xstar);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return EXIT_SUCCESS;
}

/**
 * Apply CPLEX to solve the TSP instance using Hard Fixing heuristic
 * @param inst instance
 * @param sol solution to be filled and used as starting point
 *
 * @return 0 if successful, 1 otherwise
 */
int apply_cplex_hardfix(instance *inst, solution *sol)
{
    const double start_time = second();
    double time_elapsed;
    solution current_best;
    current_best.tour = NULL;
    current_best.initialized = 0;   
    current_best.cost = CPX_INFBOUND;
    int status = EXIT_FAILURE;
    double *xstar = NULL;
    int *comp = NULL;
    static int *fixed_edges = NULL;
    static int num_fixed_edges = 0;


    // === Initialize CPLEX environment ===
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    if (!env)
    {
        print_error("CPXopenCPLEX error");
        return EXIT_FAILURE;
    }

    CPXLPptr lp = CPXcreateprob(env, &error, "TSP_HARDFIX");
    if (!lp)
    {
        print_error("CPXcreateprob error");
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // === Build initial model ===
    if (build_model(inst, env, lp))
    {
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // Store number of columns (variables)
    inst->ncols = CPXgetnumcols(env, lp);

    // === Set up callback for lazy constraints (SEC) ===
    CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
    if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst))
    {
        print_error("CPXcallbacksetfunc error");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // === Generate initial solution using greedy heuristic ===
    // Build a warm start directly
    if (sol->initialized)
    {
        copy_sol(sol, &current_best);

        if (add_warm_start(env, lp, inst, &current_best))
        {
            print_error("Failed to add warm start");
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return EXIT_FAILURE;
        }
        else if (VERBOSE >= 50)
            printf("[HARDFIX] Initial warm start cost: %.2f\n", current_best.cost);

    }
    else
    {
        // If no initial solution is provided
        solution temp_two_opt;
        temp_two_opt.tour = malloc((inst->nnodes + 1) * sizeof(int));
        if (!temp_two_opt.tour) { 
            print_error("Memory allocation failed for temp_two_opt.tour");
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return EXIT_FAILURE;
        }
        temp_two_opt.initialized = 0;

        if (apply_two_opt(inst, &temp_two_opt) != 0 || !temp_two_opt.initialized)
        {
            print_error("Two-opt warm start generation failed");
            free_sol(&temp_two_opt);
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return EXIT_FAILURE;
        }
        copy_sol(&temp_two_opt, &current_best);
        if (add_warm_start(env, lp, inst, &current_best))
        {
            print_error("Failed to add warm start");
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return EXIT_FAILURE;
        }
        else if (VERBOSE >= 50)
            printf("Warm start injected with cost %.2f\n", current_best.cost);

        free_sol(&temp_two_opt);
    }

    // === Hard fixing main loop ===
    int iteration = 0;
    int improved = 0;
    int n = inst->nnodes;

    // Allocate memory for component tracking (for potential patching)
    comp = malloc(n * sizeof(int));
    if (!comp)
    {
        print_error("Memory allocation failed for comp array");
        free_sol(&current_best);
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    if (!fixed_edges)
    {
        // Reset variable bounds to original (unfix all variables)
        fixed_edges = malloc(inst->ncols * sizeof(int));
        if (!fixed_edges)
        {
            print_error("Memory allocation failed for fixed_edges");
            free(comp);
            free_sol(&current_best);
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return EXIT_FAILURE;
        }
    }


    while (1)
    {
        iteration++;
        time_elapsed = second() - start_time;

        // Check time limit
        if (time_elapsed >= inst->timelimit)
        {
            printf("[HARDFIX] Time limit reached after %d iterations (%d improvements)\n",
                   iteration - 1, improved);
            break;
        }

        double remaining_time = inst->timelimit - time_elapsed;

        // Configure time limit for the subproblem
        CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time);
        CPXsetintparam(env, CPX_PARAM_SCRIND, VERBOSE >= 100 ? CPX_ON : CPX_OFF); // Turn on\off CPLEX output for iterations
        CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, VERBOSE >= 100 ? 4 : 0);        // Minimal display

        for (int i = 0; i < num_fixed_edges; i++) 
        {
            double lb = 0.0;
            if (CPXchgbds(env, lp, 1, &fixed_edges[i], "L", &lb))
            {
                print_error("CPXchgbds error when unfixing variables");
                goto CLEANUP;
            }
        }   
        num_fixed_edges = 0;

        // === Randomly fix edges from current_best solution ===
        int fixed_count = 0;

        // First extract the edges from the current best tour
        for (int i = 0; i < n; i++)
        {
            int from = current_best.tour[i];
            int to = current_best.tour[i + 1];

            if (from > to)
            {
                int temp = from;
                from = to;
                to = temp;
            }

            // Decide randomly whether to fix this edge
            if (((double)rand() / RAND_MAX) < inst->hard_fixing_percentage)
            {
                int idx = xpos(from, to, inst);
                double lb = 1.0; // Fix to 1 (edge must be used)

                if (CPXchgbds(env, lp, 1, &idx, "L", &lb))
                {
                    print_error("CPXchgbds error when fixing variables");
                    goto CLEANUP;
                }
                fixed_edges[num_fixed_edges++] = idx;
                fixed_count++;
            }
        }

        if (VERBOSE >= 50)
            printf("[HARDFIX] Iter %2d | Fixed edges: %2d | Remaining time: %.2fs\n", iteration, fixed_count, inst->timelimit - time_elapsed);

        // Add MIP start from current best solution
        static int current_best_changed = 1; // Set to true initially

        if (current_best_changed)
        {

            if (add_warm_start(env, lp, inst, &current_best))
            {
                fprintf(stderr, "Warning: Failed to add MIP start in iteration %d\n", iteration);
            }

            current_best_changed = 0;
        }

        // === Solve the restricted problem ===
        if (CPXmipopt(env, lp))
        {
            print_error("CPXmipopt error");
            goto CLEANUP;
        }

        // Check if we found a solution
        int sol_stat = CPXgetstat(env, lp);
        if (sol_stat == CPXMIP_OPTIMAL || sol_stat == CPXMIP_OPTIMAL_TOL ||
            sol_stat == CPXMIP_TIME_LIM_FEAS)
        {

            // Allocate or reallocate xstar if needed
            if (!xstar)
            {
                xstar = calloc(inst->ncols, sizeof(double));
                if (!xstar)
                {
                    print_error("Memory allocation failed for xstar");
                    goto CLEANUP;
                }
            }

            if (CPXgetx(env, lp, xstar, 0, inst->ncols - 1))
            {
                print_error("CPXgetx error");
                goto CLEANUP;
            }

            // Build the solution
            solution new_sol;
            new_sol.initialized = 0;
            new_sol.tour = NULL;

            if (build_solution(xstar, inst, &new_sol))
            {
                print_error("Failed to build solution");
                goto CLEANUP;
            }

            // Check if we improved
            if (new_sol.cost < current_best.cost - EPSILON)
            {
                current_best_changed = 1;

                improved++;
                if (VERBOSE >= 10)
                    printf("[HARDFIX] Iteration %d: Improved from %.2f to %.2f (-%.2f)\n",
                           iteration, current_best.cost, new_sol.cost, current_best.cost - new_sol.cost);

                free_sol(&current_best);
                copy_sol(&new_sol, &current_best);
            }
            else if (VERBOSE >= 70)
            {
                printf("[HARDFIX] Iteration %d: No improvement (current best: %.2f)\n",
                       iteration, current_best.cost);
            }

            // Adaptive tuning of fixing percentage
            if (new_sol.cost < current_best.cost - EPSILON)
            {
                inst->hard_fixing_percentage = fmin(1.0, inst->hard_fixing_percentage + 0.02);
            }
            else
            {
                inst->hard_fixing_percentage = fmax(0.1, inst->hard_fixing_percentage - 0.01);
            }


            free_sol(&new_sol);
        }
        else if (VERBOSE >= 50)
        {
            printf("[HARDFIX] Iteration %d: No feasible solution found (status %d)\n", iteration, sol_stat);
        }
    }

    // === Time expired - Check final solution and try patching if needed ===
    if (xstar)
    {
        int ncomp = build_components(xstar, comp, inst);

        // If the last solution had multiple components, try to patch it
        if (ncomp > 1)
        {
            if (VERBOSE >= 30)
                printf("[HARDFIX] Attempting to patch solution with %d components...\n", ncomp);

            solution patched_sol;
            patched_sol.initialized = 0;
            patched_sol.tour = NULL;

            if (patch_solution(xstar, inst, &patched_sol) == EXIT_SUCCESS)
            {
                // Check if patched solution is better than current best
                if (patched_sol.initialized && patched_sol.cost < current_best.cost - EPSILON)
                {
                    printf("Patched solution improved from %.2f to %.2f\n",
                           current_best.cost, patched_sol.cost);
                    free_sol(&current_best);
                    copy_sol(&patched_sol, &current_best);
                }
                free_sol(&patched_sol);
            }
        }
    }

    // Copy best solution to output
    if (current_best.initialized)
    {
        copy_sol(&current_best, sol);
        status = EXIT_SUCCESS;
    }
    else
    {
        print_error("No valid solution found during hard fixing");
        status = EXIT_FAILURE;
    }

CLEANUP:
    if (xstar)
        free(xstar);
    if (comp)
        free(comp);
    if (fixed_edges)
    {
        free(fixed_edges);
        fixed_edges = NULL;
    }

    free_sol(&current_best);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return status;
}

/**
 * Apply CPLEX to solve the TSP instance using Local Branching (Soft Fixing) heuristic
 * @param inst instance
 * @param sol solution to be filled and used as starting point
 *
 * @return 0 if successful, 1 otherwise
 */
int apply_cplex_localbranch(instance *inst, solution *sol)
{
    const double start_time = second();
    double time_elapsed;
    solution current_best;
    current_best.tour = NULL;

    int status = EXIT_FAILURE;
    double *xstar = NULL;
    int *comp = NULL;
    int n = inst->nnodes;

    // === Initialize CPLEX environment ===
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    if (!env)
    {
        print_error("CPXopenCPLEX error");
        return EXIT_FAILURE;
    }

    CPXLPptr lp = CPXcreateprob(env, &error, "TSP_LOCALBRANCH");
    if (!lp)
    {
        print_error("CPXcreateprob error");
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // === Build initial model ===
    if (build_model(inst, env, lp))
    {
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // Store number of columns (variables)
    inst->ncols = CPXgetnumcols(env, lp);

    // === Set up callback for lazy constraints (SEC) ===
    CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
    if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst))
    {
        print_error("CPXcallbacksetfunc error");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // === Generate initial solution using greedy heuristic ===
    solution initial_sol;
    initial_sol.initialized = 0;
    initial_sol.cost = CPX_INFBOUND;
    initial_sol.tour = NULL;

    initial_sol.tour = malloc((n + 1) * sizeof(int));
    if (!initial_sol.tour)
    {
        print_error("Memory allocation failed for initial_sol.tour");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    if (sol->initialized) 
    {
        copy_sol(sol, &initial_sol);
    }    
    else 
    {
        if (apply_two_opt(inst, &initial_sol) != 0 || !initial_sol.initialized)
        {
            print_error("Two opt warm start generation failed");
            free_sol(&initial_sol);
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return EXIT_FAILURE;
        }
    }

    if (add_warm_start(env, lp, inst, &initial_sol))
    {
        print_error("Failed to add warm start");
        free_sol(&initial_sol);
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }
    else if (VERBOSE >= 50)
        printf("[LOC.BRANCH] Initial warm start cost: %.2f\n", initial_sol.cost);

    // Copy as a starting point for Local Branching
    copy_sol(&initial_sol, &current_best);

    // === Local Branching main loop ===
    int k = 20;                         // Initial neighborhood size (k-opt parameter)
    const int k_min = 5;                // Minimum neighborhood size
    const int k_max = min(50, n/2);     // Maximum neighborhood size, It is used to avoid a too large neighborhood
    const int node_limit = 1000;        // Node limit per iteration
    int iteration = 0;
    int improved = 0;
    int consecutive_no_improvement = 0;

    // Allocate memory for component tracking (for potential patching)
    comp = malloc(n * sizeof(int));
    if (!comp)
    {
        print_error("Memory allocation failed for comp array");
        free_sol(&initial_sol);
        free_sol(&current_best);
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // Allocate memory for local branching constraint
    int *indices = malloc(inst->ncols * sizeof(int));
    double *values = malloc(inst->ncols * sizeof(double));
    if (!indices || !values)
    {
        print_error("Memory allocation failed for constraint arrays");
        if (indices)
            free(indices);
        if (values)
            free(values);
        free(comp);
        free_sol(&initial_sol);
        free_sol(&current_best);
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // Allocate memory for xstar
    xstar = calloc(inst->ncols, sizeof(double));
    if (!xstar)
    {
        print_error("Memory allocation failed for xstar");
        free(indices);
        free(values);
        free(comp);
        free_sol(&initial_sol);
        free_sol(&current_best);
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // Pre-allocate array to track edges in current solution
    int *edge_in_solution = calloc(inst->ncols, sizeof(int));
    if (!edge_in_solution)
    {
        print_error("Memory allocation failed for edge_in_solution");
        free(xstar);
        free(indices);
        free(values);
        free(comp);
        free_sol(&initial_sol);
        free_sol(&current_best);
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return EXIT_FAILURE;
    }

    // Constraint row number
    int lb_row = -1;

    if (VERBOSE >= 30)
        printf("[LOC.BRANCH] Starting with initial solution cost: %.2f\n", current_best.cost);

    while (1)
    {
        iteration++;
        time_elapsed = second() - start_time;

        // Check time limit
        if (time_elapsed >= inst->timelimit)
        {
            if (VERBOSE >= 30)
                printf("Time limit reached after %d iterations (%d improvements)\n",
                        iteration - 1, improved);
            break;
        }

        double remaining_time = inst->timelimit - time_elapsed;

        // Configure parameters for the subproblem
        CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time);
        CPXsetintparam(env, CPX_PARAM_NODELIM, node_limit); // Limit nodes per iteration
        CPXsetintparam(env, CPX_PARAM_SCRIND, VERBOSE >= 100 ? CPX_ON : CPX_OFF);
        CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, VERBOSE >= 100 ? 4 : 0);

        // Remove previous local branching constraint if it exists
        if (lb_row >= 0)
        {
            if (CPXdelrows(env, lp, lb_row, lb_row))
            {
                print_error("CPXdelrows error when removing local branching constraint");
                goto CLEANUP;
            }
            lb_row = -1;
        }

        // === Create local branching constraint ===
        // Reset edge tracking array
        memset(edge_in_solution, 0, inst->ncols * sizeof(int));

        // Mark edges in current solution
        for (int i = 0; i < n; i++)
        {
            int u = current_best.tour[i];
            int v = current_best.tour[i + 1];
            edge_in_solution[xpos(u, v, inst)] = 1;
        }

        // Local branching constraint (Hamming distance ≤ k):
        //    ∑_{e: x^H_e = 0} x_e  +  ∑_{e: x^H_e = 1} (1 - x_e)  ≤ k

        int nnz = 0;
        for (int e = 0; e < inst->ncols; e++)
        {
            indices[nnz] = e;
            // if edge_in_solution[e]==1 → coeff = -1 (1 - x_e)
            // otherwise coeff = +1 (x_e)
            values[nnz] = edge_in_solution[e] ? -1.0 : 1.0;
            nnz++;
        }
        double rhs   = (double)k;
        char   sense = 'L';
        int    rmatbeg[1] = { 0 };

        if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense,
                       rmatbeg, indices, values, NULL, NULL))
        {
            print_error("CPXaddrows error when adding local branching constraint");
            goto CLEANUP;
        }

        // Get the index of the newly added row
        lb_row = CPXgetnumrows(env, lp) -1;

        if (VERBOSE >= 60)
            printf("[LOC.BRANCH] Iteration %d: Added constraint with k = %d, RHS = %.0f\n", 
                   iteration, k, rhs);

        // Add MIP start from current best solution
        if (add_warm_start(env, lp, inst, &current_best))
        {
            if (VERBOSE >= 70)
                fprintf(stderr, "Warning: Failed to add MIP start in iteration %d\n", iteration);
        }

        // === Solve the restricted problem ===
        if (CPXmipopt(env, lp))
        {
            print_error("CPXmipopt error");
            goto CLEANUP;
        }

        // Check solution status
        int sol_stat = CPXgetstat(env, lp);
        int is_optimal = (sol_stat == CPXMIP_OPTIMAL || sol_stat == CPXMIP_OPTIMAL_TOL);
        int is_feasible = is_optimal || (sol_stat == CPXMIP_TIME_LIM_FEAS) || (sol_stat == CPXMIP_NODE_LIM_FEAS);
        int is_time_limit = (sol_stat == CPXMIP_TIME_LIM_FEAS || sol_stat == CPXMIP_TIME_LIM_INFEAS);

        if (VERBOSE >= 70)
            printf("[LOC.BRANCH] Iteration %d: CPLEX status = %d, feasible = %d, optimal = %d\n", 
                   iteration, sol_stat, is_feasible, is_optimal);

        if (!is_feasible)
        {
            // Problem was infeasible or no solution found
            consecutive_no_improvement++;
            
            if (VERBOSE >= 50)
                printf("[LOC.BRANCH] No feasible solution found (status = %d)\n", sol_stat);

            // Increase k when no feasible solution is found
            k = min(k + 5, k_max);
            if (VERBOSE >= 60)
                printf("[LOC.BRANCH] Increasing k to %d due to infeasibility\n", k);

            // If k reaches maximum, terminate
            if (k >= k_max)
            {
                if (VERBOSE >= 50)
                    printf("[LOC.BRANCH] Reached maximum k, terminating\n");
                break;
            }
            continue;
        }

        // Extract the solution
        if (CPXgetx(env, lp, xstar, 0, inst->ncols - 1))
        {
            print_error("CPXgetx error");
            goto CLEANUP;
        }

        // Build the solution
        solution new_sol;
        new_sol.initialized = 0;
        new_sol.tour = NULL;

        // Check if the solution is a valid tour
        int ncomp = build_components(xstar, comp, inst);
        if (ncomp > 1)
        {
            // Try to patch the solution into a valid tour
            if (patch_solution(xstar, inst, &new_sol) != EXIT_SUCCESS)
            {
                if (VERBOSE >= 60)
                    printf("[LOC.BRANCH] Patching failed for %d components\n", ncomp);
                consecutive_no_improvement++;
                
                // Adjust k and continue
                if (is_optimal)
                    k = min(k + 3, k_max);
                else if (is_time_limit)
                    k = max(k - 2, k_min);
                
                continue;
            }

            if (VERBOSE >= 50)
                printf("[LOC.BRANCH] Patched disconnected solution (%d components)\n", ncomp);
        }
        else
        {
            // Build solution from valid tour
            if (build_solution(xstar, inst, &new_sol) != EXIT_SUCCESS)
            {
                print_error("Failed to build solution");
                consecutive_no_improvement++;
                continue;
            }
        }

        // Check if we improved
        if (new_sol.cost < current_best.cost - EPSILON)
        {
            improved++;
            consecutive_no_improvement = 0;
            
            if (VERBOSE >= 20)
                printf("[LOC.BRANCH] Iter %d: Improved from %.2f to %.2f (k=%d)\n", 
                       iteration, current_best.cost, new_sol.cost, k);

            free_sol(&current_best);
            copy_sol(&new_sol, &current_best);

            // Reset k after improvement
            k = 20; // Reset to initial value, not k_min
        }
        else
        {
            consecutive_no_improvement++;
            
            if (VERBOSE >= 70)
                printf("[LOC.BRANCH] Iteration %d: No improvement (current: %.2f, found: %.2f)\n",
                       iteration, current_best.cost, new_sol.cost);

            // Adjust k based on the solution status
            if (is_optimal)
            {
                // Solution is optimal but no improvement -> increase k
                k = min(k + 3, k_max);
                if (VERBOSE >= 60)
                    printf("[LOC.BRANCH] Optimal without improvement, increasing k to %d\n", k);
            }
            else if (is_time_limit)
            {
                // Time limit reached -> decrease k
                k = max(k - 2, k_min);
                if (VERBOSE >= 60)
                    printf("[LOC.BRANCH] Time limit reached, decreasing k to %d\n", k);
            }
        }

        free_sol(&new_sol);

        // Termination conditions
        if (k >= k_max && consecutive_no_improvement >= 3)
        {
            if (VERBOSE >= 30)
                printf("[LOC.BRANCH] Reached maximum k with no improvement, terminating\n");
            break;
        }

        if (consecutive_no_improvement >= 10)
        {
            if (VERBOSE >= 30)
                printf("[LOC.BRANCH] Too many iterations without improvement, terminating\n");
            break;
        }
    }

    // Copy best solution to output
    copy_sol(&current_best, sol);
    status = EXIT_SUCCESS;

    if (VERBOSE >= 20)
        printf("[LOC.BRANCH] Completed: %d iterations, %d improvements, final cost: %.2f\n", 
               iteration - 1, improved, current_best.cost);

CLEANUP:
    if (xstar) free(xstar);
    if (comp) free(comp);
    if (indices) free(indices);
    if (values) free(values);
    if (edge_in_solution) free(edge_in_solution);
    free_sol(&current_best);
    free_sol(&initial_sol);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return status;
}
