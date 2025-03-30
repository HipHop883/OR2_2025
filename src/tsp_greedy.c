#include "utils.h"
#include "tsp.h"
#include "tsp_greedy.h"

#define N_GREEDY_STARTS 10 // Number of random starts for multi-start greedy

/**
 * Applies a multi-start greedy heuristic using the nearest neighbor method.
 * Tries N_GREEDY_STARTS random starting points and keeps the best result.
 * Stops early if time limit is reached.
 *
 * @param inst pointer to the TSP instance
 * @param best_sol solution structure to be filled with the best found tour
 * @return 0 on success
 */
int apply_greedy_search(const instance *inst, solution *best_sol)
{
    set_seed(inst->randomseed);

    solution current_sol;
    current_sol.tour = (int *)malloc((inst->nnodes + 1) * sizeof(int));

    if (!current_sol.tour)
    {
        print_error("Memory allocation failed in solve_greedy");
        return EXIT_FAILURE;
    }

    double best_cost = CPX_INFBOUND;

    if (best_sol->initialized)      // Check if already initialized
    {
        best_cost = best_sol->cost; // Use the existing cost
    }

    for (int i = 0; i < N_GREEDY_STARTS; i++)
    {
        if (check_time(inst))
        {
            if (VERBOSE >= 20)
                printf("Greedy search stopped early due to time limit.\n");
            break;
        }

        int start_node = rand() % inst->nnodes;
        apply_nearest_neighbor(inst, &current_sol, start_node);

        if (current_sol.cost < best_cost)
        {
            best_cost = current_sol.cost;
            memcpy(best_sol->tour, current_sol.tour, (inst->nnodes + 1) * sizeof(int));
        }
    }

    best_sol->cost = best_cost;   // Update best_sol cost
    best_sol->initialized = 1;    // Mark as initialized

    free(current_sol.tour);
    return EXIT_SUCCESS;
}

/**
 * Builds a TSP tour starting from a given node using the nearest neighbor heuristic.
 *
 * @param inst pointer to the TSP instance
 * @param sol solution structure to fill with the constructed tour
 * @param start_node index of the starting node
 * @return 0 on success
 */
int apply_nearest_neighbor(const instance *inst, solution *sol, int start_node)
{
    double time = second();
    if (VERBOSE >= 50)
    {
        printf("Applying nearest neighbor heuristic from node %d...\n", start_node);
    }

    int nnodes = inst->nnodes;
    int *visited = (int *)calloc(nnodes, sizeof(int));
    int current = start_node;
    visited[current] = 1;
    sol->tour[0] = current;

    for (int i = 1; i < nnodes; i++)
    {
        double min_cost = CPX_INFBOUND;
        int next = -1;

        for (int j = 0; j < nnodes; j++)
        {
            if (!visited[j])
            {
                double c = inst->cost_matrix[current][j];
                if (c < min_cost)
                {
                    min_cost = c;
                    next = j;
                }
            }
        }

        visited[next] = 1;
        sol->tour[i] = next;
        current = next;
    }

    sol->tour[nnodes] = sol->tour[0]; // close the tour
    free(visited);

    if (VERBOSE >= 50)
    {
        printf("Nearest neighbor heuristic applied\n");
        printf("Elapsed time: %lf\n", second() - time);
        printf("--------------------------------------------\n");
    }

    evaluate_path_cost(inst, sol);

    if (!sol->initialized)      // Mark as initialized only if not already initialized
    {
        sol->initialized = 1;
    }

    return EXIT_SUCCESS;
}
