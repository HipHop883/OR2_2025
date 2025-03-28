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
        return 1;

    solution *current_sol = (solution *)malloc(sizeof(solution));
    solution *best_sol = (solution *)malloc(sizeof(solution));

    current_sol->tour = (int *)malloc(sizeof(int) * (tsp->nnodes + 1));
    best_sol->tour = (int *)malloc(sizeof(int) * (tsp->nnodes + 1));

    if (!current_sol->tour || !best_sol->tour)
    {
        print_error("Memory allocation failed");
        return EXIT_FAILURE;
    }

    current_sol->cost = LARGE;
    best_sol->cost = LARGE;

    // tsp->starting_time = second();

    // Get an initial solution using the greedy algorithm.
    // we run 2opt in the while loop.
    if (nearest_neighbor(tsp, current_sol, 0) != 0)
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

        cost_path(tsp, current_sol);
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
 *
 * @param inst instance
 * @param sol solution
 * @return 0 if the solution is found, 1 otherwise
 */
int tsp_solve_tabu(instance *inst, solution *sol)
{
    if (!inst || !sol)
    {
        return 1;
    }

    int nnodes = inst->nnodes;
    tabuList *tabu = init_tabu_list(nnodes);
    if (!tabu)
        return 1;

    if (nearest_neighbor(inst, sol, 0) != 0 || two_opt(inst, sol) != 0)
    {
        print_error("Nearest neighbor heuristic failed in tabu");
        free_tabu(tabu);
        return 1;
    }

    // Copy the initial solution to the best solution
    solution *best_sol = (solution *)malloc(sizeof(solution));
    best_sol->tour = (int *)malloc((nnodes + 1) * sizeof(int));
    if (!best_sol->tour || !best_sol)
    {
        print_error("Memory allocation failed");
        free_tabu(tabu);
        free(best_sol);
        return 1;
    }
    memcpy(best_sol->tour, sol->tour, (nnodes + 1) * sizeof(int));
    best_sol->cost = sol->cost;

    // Copy the initial solution to the current solution
    solution *current_sol = (solution *)malloc(sizeof(solution));
    current_sol->tour = (int *)malloc((nnodes + 1) * sizeof(int));
    if (!current_sol->tour || !current_sol)
    {
        print_error("Memory allocation failed");
        free_tabu(tabu);
        free(best_sol->tour);
        free(best_sol);
        free(current_sol);
        return 1;
    }
    memcpy(current_sol->tour, sol->tour, (nnodes + 1) * sizeof(int));
    current_sol->cost = sol->cost;
    printf("Inital cost: %lf\n", best_sol->cost);

    // Initialize the number of iterations
    int iter = 0;

    while (second() - inst->starting_time < inst->timelimit)
    {
        if (best_2opt_not_tabu(inst, current_sol, tabu) == -1)
            print_error("No valid 2-opt move found.");

        if (current_sol->cost < best_sol->cost)
        {
            memcpy(best_sol->tour, current_sol->tour, (nnodes + 1) * sizeof(int));
            best_sol->cost = current_sol->cost;
            printf("[UPDATE] New best at iter %d: %.2lf\n", iter, best_sol->cost);
        }

        iter++;
    }

    // Copy the best found solution to the output parameters.
    memcpy(sol->tour, best_sol->tour, sizeof(int) * (inst->nnodes + 1));
    sol->cost = best_sol->cost;

    // free memory
    free_tabu(tabu);
    free(best_sol->tour);
    free(best_sol);
    free(current_sol->tour);
    free(current_sol);

    return 0;
}

/**
 * Best 2-opt move that is not tabu.
 */
static int best_2opt_not_tabu(instance *tsp, solution *sol, tabuList *tabu)
{
    int nnodes = tsp->nnodes;
    double best_delta = LARGE;
    int best_i = -1, best_j = -1;

    for (int i = 1; i < nnodes - 1; i++)
    {
        for (int j = i + 1; j < nnodes; j++)
        {
            double actual_delta = delta(i, j, sol, tsp);

            // No need to normalize here — it's handled inside is_tabu_move()
            int is_tabu = is_tabu_move(tabu, i, j);
            double new_cost = sol->cost + actual_delta;

            if (!is_tabu && actual_delta != 0 && actual_delta < best_delta)
            {
                best_delta = actual_delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if (best_i == -1 || best_j == -1)
        return -1;

    swap_path(best_i, best_j, sol);

    if (cost_path(tsp, sol))
        return -1;

    add_tabu_move(tabu, best_i, best_j);
    return 0;
}

/**
 * Check if move is in tabu list.
 */
int is_tabu_move(tabuList *tabu, int i, int j)
{
    if (i > j)
    {
        int tmp = i;
        i = j;
        j = tmp;
    }
    for (int k = 0; k < tabu->size; k++)
        if (tabu->tabu_list[k][0] == i && tabu->tabu_list[k][1] == j)
            return 1;
    return 0;
}

/**
 * Add move to tabu list (FIFO).
 */
void add_tabu_move(tabuList *tabu, int i, int j)
{
    if (i > j)
    {
        int tmp = i;
        i = j;
        j = tmp;
    }

    if (tabu->size < tabu->tenure)
    {
        tabu->tabu_list[tabu->size][0] = i;
        tabu->tabu_list[tabu->size][1] = j;
        tabu->size++;
    }
    else
    {
        free(tabu->tabu_list[0]);
        for (int k = 1; k < tabu->tenure; k++)
            tabu->tabu_list[k - 1] = tabu->tabu_list[k];

        tabu->tabu_list[tabu->tenure - 1] = (int *)calloc(2, sizeof(int));
        tabu->tabu_list[tabu->tenure - 1][0] = i;
        tabu->tabu_list[tabu->tenure - 1][1] = j;
    }
}

/**
 * Free all memory associated with the tabu list.
 */
void free_tabu(tabuList *tabu)
{
    for (int i = 0; i < tabu->tenure; i++)
        free(tabu->tabu_list[i]);
    free(tabu->tabu_list);
    free(tabu);
}

tabuList *init_tabu_list(int nnodes)
{
    tabuList *tabu = (tabuList *)malloc(sizeof(tabuList));
    if (!tabu)
    {
        print_error("Failed to allocate tabuList");
        return NULL;
    }

    tabu->tenure = nnodes + nnodes / 100;
    tabu->tabu_list = (int **)calloc(tabu->tenure, sizeof(int *));
    if (!tabu->tabu_list)
    {
        print_error("Failed to allocate tabu list");
        free(tabu);
        return NULL;
    }

    for (int i = 0; i < tabu->tenure; i++)
    {
        tabu->tabu_list[i] = (int *)calloc(2, sizeof(int));
        if (!tabu->tabu_list[i])
        {
            print_error("Failed to allocate tabu entry");
            for (int j = 0; j < i; j++)
                free(tabu->tabu_list[j]);
            free(tabu->tabu_list);
            free(tabu);
            return NULL;
        }
    }

    tabu->size = 0;
    return tabu;
}
