#include "utils.h"
#include "tsp.h"
#include "tsp_greedy.h"

#define DIVIDER "-----------------------------------------------------------"

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
    double starting_time = second();

    solution current_sol;
    current_sol.tour = NULL;
    current_sol.cost = 0.0;
    current_sol.initialized = 0;
    
    current_sol.tour = (int *)malloc((inst->nnodes + 1) * sizeof(int));
    if (!current_sol.tour)
    {
        print_error("Memory allocation failed in solve_greedy");
        return EXIT_FAILURE;
    }

    double best_cost = CPX_INFBOUND;
    if (best_sol->initialized)
    {
        best_cost = best_sol->cost;
    }

    if (!best_sol->tour)
    {
        best_sol->tour = malloc((inst->nnodes + 1) * sizeof(int));
        if (!best_sol->tour)
        {
            print_error("Memory allocation failed for best_sol->tour");
            free(current_sol.tour);
            return EXIT_FAILURE;
        }
    }

    int greedy_starts = fmin(inst->greedy_starts, inst->nnodes);

    for (int i = 0; i < greedy_starts; i++)
    {
        if (check_time(inst, starting_time))
        {
            if (VERBOSE >= 10)
                printf("[GREEDY] %-38s%11d\n", "Stopped early due to time limit:", i);
            break;
        }

        int start_node = rand() % inst->nnodes;

        if (VERBOSE >= 50)
        {
            printf("[GREEDY] %-38s%11d\n", "Start index:", i + 1);
            printf("[GREEDY] %-38s%11d\n", "Starting node:", start_node);
        }

        apply_nearest_neighbor(inst, &current_sol, start_node);

        if (VERBOSE >= 20)
            printf("[GREEDY] %-38s%11.2f\n", "Tour cost from node:", current_sol.cost);

        if (current_sol.cost < best_cost)
        {
            best_cost = current_sol.cost;
            memcpy(best_sol->tour, current_sol.tour, (inst->nnodes + 1) * sizeof(int));

            if (VERBOSE >= 10)
                printf("[GREEDY] %-38s%11.2f\n", "NEW BEST FOUND. Cost:", best_cost);
        }

        if (VERBOSE >= 50)
            printf("[GREEDY] %-38s%11.2f\n", "Current best after this start:", best_cost);

        printf("%s\n", DIVIDER);
    }

    best_sol->cost = best_cost;
    best_sol->initialized = 1;

    if (VERBOSE >= 10)
    {
        printf("[GREEDY] %-38s%11.2f\n", "Final best cost:", best_cost);
        printf("[GREEDY] %-38s%11.3fs\n", "Total elapsed time:", second() - starting_time);
    }

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

    int nnodes = inst->nnodes;
    int *visited = (int *) calloc(nnodes, sizeof(int));
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
        if (next == -1) 
        {
            fprintf(stderr, "[ERROR] Nearest neighbor: no unvisited nodes found at i = %d\n", i);
            free(visited);
            return EXIT_FAILURE;
        }
        visited[next] = 1;
        sol->tour[i] = next;
        current = next;
    }

    sol->tour[nnodes] = sol->tour[0];
    free(visited);

    evaluate_path_cost(inst, sol);

    if (VERBOSE >= 50)
        printf("[NN]     %-38s%11.3fs\n", "Tour built in:", second() - time);

    if (!sol->initialized)
        sol->initialized = 1;

    return EXIT_SUCCESS;
}
