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

/**
 * Generate 3-opt positions: fills positions[0..2] with valid indices for a 3-opt move.
 * @param tsp TSP instance
 * @param positions array of 3 integers
 * @return void
 */
static void generate_3opt_positions(instance *tsp, int *positions)
{
    while (1)
    {
        positions[0] = rand() % tsp->nnodes;
        positions[1] = rand() % tsp->nnodes;
        positions[2] = rand() % tsp->nnodes;
        qsort(positions, 3, sizeof(int), compar);

        // Ensure the positions are sufficiently apart:
        if (positions[1] > positions[0] + 1 &&
            positions[2] > positions[1] + 1 &&
            positions[2] + 1 < tsp->nnodes)
        {
            break;
        }
    }
}

/**
 * Compare two integers for qsort
 * @param a integer a
 * @param b integer b
 * @return the difference between a and b
 */
static int compar(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

/**
 * Recompute the solution cost after applying a 3-opt swap.
 * @param tsp TSP instance
 * @param solution solution path
 * @return 0 if the solution cost is recomputed successfully, -1 otherwise
 */
int tsp_3opt_solution(instance *tsp, solution *sol)
{
    if (!tsp->cost_matrix || tsp->nnodes <= 0)
    {
        print_error("Invalid TSP instance in tsp_3opt_solution");
        return -1;
    }
    if (VERBOSE >= 50)
    {
        printf("Applying three-opt solution...\n");
    }
    int positions[3];
    generate_3opt_positions(tsp, positions);

    int *temp_tour = (int *)malloc(sizeof(int) * (tsp->nnodes + 1));
    if (temp_tour == NULL)
    {
        print_error("Memory allocation failed in tsp_3opt_solution");
        return -1;
    }

    // Apply the 3-opt swap using the generated positions.
    tsp_3opt_swap(tsp, sol, temp_tour, positions);
    memcpy(sol->tour, temp_tour, sizeof(int) * (tsp->nnodes + 1));

    free(temp_tour);

    // Recompute and update the solution cost.
    sol->cost = cost_path(tsp, sol);

    printf("Three-opt solution applied\n");
    printf("--------------------------------------------\n");
    return 0;
}

/**
 * Perform a 3-opt swap on the current solution, storing the result in the new solution.
 * @param tsp TSP instance
 * @param current_sol current solution
 * @param new_tour new solution path
 * @param positions array of 3 integers
 * @return void
 */
static void tsp_3opt_swap(instance *tsp, solution *current_sol, int *new_tour, int *positions)
{
    // Copy segment from index 0 to positions[0] unchanged.
    memcpy(new_tour, current_sol->tour, sizeof(int) * (positions[0] + 1));
    int pos = positions[0] + 1;

    // Reverse the segment between positions[0]+1 and positions[1].
    for (int i = positions[1]; i >= positions[0] + 1; i--)
    {
        new_tour[pos++] = current_sol->tour[i];
    }

    // Reverse the segment between positions[1]+1 and positions[2].
    for (int i = positions[2]; i >= positions[1] + 1; i--)
    {
        new_tour[pos++] = current_sol->tour[i];
    }

    // Copy the remainder of the tour from positions[2]+1 to the end.
    memcpy(new_tour + pos, current_sol->tour + positions[2] + 1, sizeof(int) * (tsp->nnodes - positions[2] - 1));
}
