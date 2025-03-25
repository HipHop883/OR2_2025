#include "utils.h"
#include "tsp.h"
#include "tsp_greedy.h"

int solve_greedy(const instance *inst, solution *best_sol)
{
    solution current_sol;
    current_sol.tour = (int *)malloc((inst->nnodes + 1) * sizeof(int));

    double best_cost = CPX_INFBOUND;

    for (int i = 0; i < inst->nnodes; i++)
    {
        int start_node = rand() % inst->nnodes;
        nearest_neighbor(inst, &current_sol, start_node);

        cost_path(inst, &current_sol);
        if (current_sol.cost < best_cost)
        {
            best_cost = current_sol.cost;
            memcpy(best_sol->tour, current_sol.tour, (inst->nnodes + 1) * sizeof(int));
        }
    }

    free(current_sol.tour);

    return EXIT_SUCCESS;
}

/**
 * Nearest neighbor heuristic
 * @param inst instance
 * @param sol solution path
 * @return 0 if the nearest neighbor heuristic is applied successfully, 1 otherwise
 */
int nearest_neighbor(const instance *inst, solution *sol, int start_node)
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
        { // find the nearest neighbor
            if (visited[j] == 0)
            {
                double c = inst->cost_matrix[current][j];
                if (c < min_cost)
                {                 // update the nearest neighbor
                    min_cost = c; // update the minimum cost
                    next = j;     // update the next node
                }
            }
        }
        visited[next] = 1;
        sol->tour[i] = next;
        current = next;
    }
    sol->tour[nnodes] = sol->tour[0]; // Close the tour by returning to the starting node
    free(visited);

    if (VERBOSE >= 50)
    {
        printf("Nearest neighbor heuristic applied\n");
        printf("Elapsed time: %lf\n", second() - time);
        printf("--------------------------------------------\n");
    }

    cost_path(inst, sol);

    return EXIT_SUCCESS;
}
