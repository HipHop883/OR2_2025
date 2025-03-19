#include <../include/vns.h>

/**
 * Solve the TSP instance using the VNS algorithm. Starting from a greedy solution, we apply
 * intensification (2-opt improvement) and diversification (3-opt kicks) until the time limit is reached.
 * @param tsp TSP instance
 * @param sol solution
 * @return 0 if the solution is found, 1 otherwise
 */
int tsp_solve_vns(const instance *tsp, solution *sol)
{
    if (!tsp || !sol->tour || !sol->cost  || tsp->nnodes <= 0 || !tsp->cost_matrix) 
    {
        return 1;
    }

    int *current_solution = (int *)malloc(sizeof(int) * tsp->nnodes);
    int *best_solution = (int *)malloc(sizeof(int) * tsp->nnodes);

    if (!current_solution || !best_solution)
    {
        print_error("Memory allocation failed");
        free(current_solution);
        current_solution = NULL;
        free(best_solution);
        best_solution = NULL;
        return 1;
    }

    double current_cost = 1e30;
    double best_cost = 1e30;

    //tsp->starting_time = second();

    // Get an initial solution using the greedy algorithm.
    // we run 2opt in the while loop.
    if (tsp_solve_greedy(tsp, current_solution, &current_cost, 0) != 0)
    {
        free(current_solution);
        free(best_solution);
        return -1;
    }

    memcpy(best_solution, current_solution, sizeof(int) * tsp->nnodes);
    best_cost = current_cost;

    //tsp->starting_time = second();

    // VNS parameters for diversification: number of 3-opt kicks.
    int min_kicks = 1;
    int max_kicks = 5;

    while (1)
    {
        if (tsp->timelimit > 0 && (tsp->starting_time + tsp->timelimit) < second())
        {
            break;
        }

        // Intensification phase: apply 2-opt improvement until no further improvement is found.
        if (tsp_2opt_solution(tsp, current_solution, &current_cost) == 1)
        {
            break;
        }

        // Update the best solution if the current solution is improved.
        if (current_cost < best_cost)
        {
            best_cost = current_cost;
            memcpy(best_solution, current_solution, sizeof(int) * tsp->nnodes);
        }

        // Diversification phase: apply a random number of 3-opt kicks.
        int kicks = (rand() % (max_kicks - min_kicks + 1)) + min_kicks;
        for (int i = 0; i < kicks; i++)
        {
            if (tsp_3opt_solution(tsp, current_solution, &current_cost) != 0)
            {
                break;
            }
        }

        current_cost = tsp_recompute_solution_arg(tsp, current_solution);
    }

    // Copy the best found solution to the output parameters.
    memcpy(sol->tour, best_solution, sizeof(int) * tsp->nnodes);
    sol->cost = best_cost;

    free(current_solution);
    current_solution = NULL;
    free(best_solution);
    best_solution = NULL;

    return 0;
}