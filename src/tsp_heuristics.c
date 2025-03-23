#include "tsp_heuristics.h"
#include "utils.h"

#define LARGE 1e30

/**
 * Solve the TSP instance using the VNS algorithm. Starting from a greedy solution, we apply
 * intensification (2-opt improvement) and diversification (3-opt kicks) until the time limit is reached.
 * @param tsp TSP instance
 * @param sol solution
 * @return 0 if the solution is found, 1 otherwise
 */
int tsp_solve_vns(instance *tsp, solution *sol)
{
    if (!tsp || !sol->tour || tsp->nnodes <= 0 || !tsp->cost_matrix)
    {
        return 1;
    }

    solution *current_sol = (solution *)malloc(sizeof(solution));
    solution *best_sol = (solution *)malloc(sizeof(solution));

    current_sol->tour = (int *)malloc(sizeof(int) * (tsp->nnodes + 1));
    best_sol->tour = (int *)malloc(sizeof(int) * (tsp->nnodes + 1));

    if (current_sol->tour == NULL || best_sol->tour == NULL)
    {
        print_error("Memory allocation failed");
        return EXIT_FAILURE;
    }

    current_sol->cost = LARGE;
    best_sol->cost = LARGE;

    // tsp->starting_time = second();

    // Get an initial solution using the greedy algorithm.
    // we run 2opt in the while loop.
    if (solve_greedy(tsp, current_sol) != 0)
    {
        print_error("Nearest neighbor heuristic failed in vns");
        free(current_sol->tour);
        free(current_sol);
        free(best_sol->tour);
        free(best_sol);
        return EXIT_FAILURE;
    }

    memcpy(best_sol->tour, current_sol->tour, sizeof(int) * (tsp->nnodes + 1));
    best_sol->cost = current_sol->cost;

    // VNS parameters for diversification: number of 3-opt kicks.
    int min_kicks = 1;
    int max_kicks = 5;

    while (1)
    {
        if (tsp->timelimit > 0 && (tsp->starting_time + tsp->timelimit) < second())
        {
            printf("Time limit reached vns\n");
            break;
        }

        // Intensification phase: apply 2-opt improvement until no further improvement is found.
        if (two_opt(tsp, current_sol) != 0)
        {
            print_error("2-opt failed in vns");
        }

        // Update the best solution if the current solution is improved.
        if (current_sol->cost < best_sol->cost)
        {
            best_sol->cost = current_sol->cost;
            memcpy(best_sol->tour, current_sol->tour, sizeof(int) * (tsp->nnodes + 1));
        }

        // Diversification phase: apply a random number of 3-opt kicks.
        int kicks = (rand() % (max_kicks - min_kicks + 1)) + min_kicks;
        for (int i = 0; i < kicks; i++)
        {
            if (tsp_3opt_solution(tsp, current_sol) != 0)
            {
                print_error("3-opt failed in vns");
            }
        }

        current_sol->cost = cost_path(tsp, current_sol);
    }

    // Copy the best found solution to the output parameters.
    memcpy(sol->tour, best_sol->tour, sizeof(int) * (tsp->nnodes + 1));
    sol->cost = best_sol->cost;

    // Free memory
    free(current_sol->tour);
    free(current_sol);
    free(best_sol->tour);
    free(best_sol);

    return 0;
}
