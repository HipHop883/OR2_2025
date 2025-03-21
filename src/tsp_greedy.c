#include "chrono.h"
#include "tsp.h"
#include "tsp_greedy.h"

/**
 * Nearest neighbor heuristic
 * @param inst instance
 * @param sol solution path
 * @return 0 if the nearest neighbor heuristic is applied successfully, 1 otherwise
 */
int nearest_neighbor(const instance *inst, solution *sol)
{
    double time = second();
    if (VERBOSE >= 50)
    {
        printf("Applying nearest neighbor heuristic...\n");
    }

    int nnodes = inst->nnodes;                         // number of nodes
    int *visited = (int *)calloc(nnodes, sizeof(int)); // visited nodes
    int current = 0;                                   // current node
    visited[current] = 1;                              // mark the current node as visited
    sol->tour[0] = current;                            // start from the current node

    for (int i = 1; i < nnodes; i++)
    {
        double min_cost = CPX_INFBOUND;
        int next = -1;
        for (int j = 0; j < nnodes; j++)
        { // find the nearest neighbor
            if (visited[j] == 0)
            {
                // double c = inst->cost_matrix[flatten_coords(current, j, nnodes)];
                double c = inst->cost_matrix[current][j];
                if (c < min_cost)
                {                 // update the nearest neighbor
                    min_cost = c; // update the minimum cost
                    next = j;     // update the next node
                }
            }
        }
        visited[next] = 1;   // mark the next node as visited
        sol->tour[i] = next; // update the best solution
        current = next;      // update the current node
    }
    sol->tour[nnodes] = sol->tour[0]; // return to the initial node
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
