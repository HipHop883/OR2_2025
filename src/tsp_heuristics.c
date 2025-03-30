#include "tsp_heuristics.h"
#include "utils.h"

static int is_tabu_move(tabuList *tabu, int i, int j)
{
    if (i > j)
    {
        int t = i;
        i = j;
        j = t;
    }

    for (int k = 0; k < tabu->size; k++)
        if (tabu->tabu_list[k][0] == i && tabu->tabu_list[k][1] == j)
            return 1;

    return 0;
}

static void add_tabu_move(tabuList *tabu, int i, int j)
{
    if (i > j)
    {
        int t = i;
        i = j;
        j = t;
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

static void free_tabu(tabuList *tabu)
{
    for (int i = 0; i < tabu->tenure; i++)
        free(tabu->tabu_list[i]);
    free(tabu->tabu_list);
    free(tabu);
}

static int best_2opt_not_tabu(instance *inst, solution *sol, tabuList *tabu)
{
    int nnodes = inst->nnodes;
    double best_delta = CPX_INFBOUND;
    int best_i = -1, best_j = -1;

    for (int i = 1; i < nnodes - 1; i++)
    {
        for (int j = i + 1; j < nnodes; j++)
        {
            double delta = path_cost_delta(i, j, sol, inst);
            int is_tabu = is_tabu_move(tabu, i, j);
            double new_cost = sol->cost + delta;

            if (!is_tabu && delta < best_delta && delta != 0)
            {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if (best_i == -1 || best_j == -1)
        return -1;

    reverse_path_segment(best_i, best_j, sol);

    if (evaluate_path_cost(inst, sol) != 0)
        return -1;

    add_tabu_move(tabu, best_i, best_j);
    return 0;
}

static tabuList *init_tabu_list(int nnodes)
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
            for (int j = 0; j < i; j++)
                free(tabu->tabu_list[j]);
            free(tabu->tabu_list);
            free(tabu);
            print_error("Failed to allocate tabu entry");
            return NULL;
        }
    }

    tabu->size = 0;
    return tabu;
}

/**
 * Solve the TSP instance using the VNS algorithm. Starting from a greedy solution, we apply
 * intensification (2-opt improvement) and diversification (3-opt kicks) until the time limit is reached.
 * @param tsp TSP instance
 * @param sol solution
 * @return 0 if the solution is found, 1 otherwise
 */
int apply_heuristic_vns(instance *inst, solution *sol)
{
    if (!inst || !sol->tour || inst->nnodes <= 0 || !inst->cost_matrix)
        return EXIT_FAILURE;

    solution *current_sol = (solution *)malloc(sizeof(solution));
    solution *best_sol = (solution *)malloc(sizeof(solution));

    current_sol->tour = (int *)malloc(sizeof(int) * (inst->nnodes + 1));
    best_sol->tour = (int *)malloc(sizeof(int) * (inst->nnodes + 1));

    if (!current_sol->tour || !best_sol->tour)
    {
        print_error("Memory allocation failed");
        return EXIT_FAILURE;
    }

    current_sol->cost = CPX_INFBOUND;
    best_sol->cost = CPX_INFBOUND;

    if (apply_nearest_neighbor(inst, current_sol, 0) != 0)
    {
        print_error("Nearest neighbor heuristic failed in VNS");
        free(current_sol->tour);
        free(current_sol);
        free(best_sol->tour);
        free(best_sol);
        return EXIT_FAILURE;
    }

    memcpy(best_sol->tour, current_sol->tour, sizeof(int) * (inst->nnodes + 1));
    best_sol->cost = current_sol->cost;

    int min_kicks = inst->vns_kmin;
    int max_kicks = inst->vns_kmax;

    while (!check_time(inst))
    {
        if (apply_two_opt(inst, current_sol) != 0)
            print_error("2-opt failed in VNS");

        if (current_sol->cost < best_sol->cost)
        {
            best_sol->cost = current_sol->cost;
            memcpy(best_sol->tour, current_sol->tour, sizeof(int) * (inst->nnodes + 1));
        }

        int kicks = (rand() % (max_kicks - min_kicks + 1)) + min_kicks;
        for (int i = 0; i < kicks; i++)
        {
            if (apply_three_opt(inst, current_sol) != 0)
                print_error("3-opt failed in VNS");
        }

        evaluate_path_cost(inst, current_sol);
    }

    if (VERBOSE >= 20)
        printf("VNS stopped due to time limit.\n");

    memcpy(sol->tour, best_sol->tour, sizeof(int) * (inst->nnodes + 1));
    sol->cost = best_sol->cost;

    free(current_sol->tour);
    free(current_sol);
    free(best_sol->tour);
    free(best_sol);

    return EXIT_SUCCESS;
}

int apply_heuristic_tabu(instance *inst, solution *sol)
{
    if (!inst || !sol)
        return EXIT_FAILURE;

    int nnodes = inst->nnodes;
    tabuList *tabu = init_tabu_list(nnodes);
    if (!tabu)
        return EXIT_FAILURE;

    if (apply_nearest_neighbor(inst, sol, 0) != 0 || apply_two_opt(inst, sol) != 0)
    {
        print_error("Initial construction failed in Tabu Search");
        free_tabu(tabu);
        return EXIT_FAILURE;
    }

    solution *best_sol = (solution *)malloc(sizeof(solution));
    solution *current_sol = (solution *)malloc(sizeof(solution));
    best_sol->tour = (int *)malloc((nnodes + 1) * sizeof(int));
    current_sol->tour = (int *)malloc((nnodes + 1) * sizeof(int));

    if (!best_sol || !current_sol || !best_sol->tour || !current_sol->tour)
    {
        print_error("Memory allocation failed");
        free_tabu(tabu);
        if (best_sol)
        {
            free(best_sol->tour);
            free(best_sol);
        }
        if (current_sol)
        {
            free(current_sol->tour);
            free(current_sol);
        }
        return EXIT_FAILURE;
    }

    memcpy(best_sol->tour, sol->tour, sizeof(int) * (nnodes + 1));
    memcpy(current_sol->tour, sol->tour, sizeof(int) * (nnodes + 1));
    best_sol->cost = current_sol->cost = sol->cost;

    int iter = 0;

    while (!check_time(inst))
    {
        if (best_2opt_not_tabu(inst, current_sol, tabu) == -1)
            print_error("No valid 2-opt move found");

        if (current_sol->cost < best_sol->cost)
        {
            memcpy(best_sol->tour, current_sol->tour, sizeof(int) * (nnodes + 1));
            best_sol->cost = current_sol->cost;

            if (VERBOSE >= 20)
                printf("[UPDATE] Iter %d — new best cost: %.2lf\n", iter, best_sol->cost);
        }

        iter++;
    }

    if (VERBOSE >= 20)
        printf("Tabu Search stopped due to time limit.\n");

    memcpy(sol->tour, best_sol->tour, sizeof(int) * (nnodes + 1));
    sol->cost = best_sol->cost;

    // free_tabu(tabu);
    free(best_sol->tour);
    free(best_sol);
    free(current_sol->tour);
    free(current_sol);

    return EXIT_SUCCESS;
}
