#include "tsp_heuristics.h"
#include "tsp_greedy.h"
#include "utils.h"

#define DIVIDER "-----------------------------------------------------------"

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
            return EXIT_FAILURE;
    return EXIT_SUCCESS;
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

static int resize_tabu_list(tabuList *tabu, int new_tenure)
{
    if (new_tenure == tabu->tenure)
        return 0;

    // Shrink: free extra entries
    if (new_tenure < tabu->tenure)
    {
        for (int i = new_tenure; i < tabu->tenure; i++)
        {
            free(tabu->tabu_list[i]);
        }
    }

    // Resize the array
    int **new_list = realloc(tabu->tabu_list, sizeof(int *) * new_tenure);
    if (!new_list)
    {
        print_error("Failed to realloc tabu list");
        return -1;
    }
    tabu->tabu_list = new_list;

    // Expand: allocate new entries
    if (new_tenure > tabu->tenure)
    {
        for (int i = tabu->tenure; i < new_tenure; i++)
        {
            tabu->tabu_list[i] = (int *)calloc(2, sizeof(int));
            if (!tabu->tabu_list[i])
            {
                print_error("Failed to allocate new tabu entry");
                return -1;
            }
        }
    }

    // Clamp current size
    if (tabu->size > new_tenure)
        tabu->size = new_tenure;

    tabu->tenure = new_tenure;
    return 0;
}

static void free_tabu(tabuList *tabu)
{
    if (!tabu)
        return;

    if (tabu->tabu_list)
    {
        for (int i = 0; i < tabu->tenure; i++)
        {
            if (tabu->tabu_list[i])
                free(tabu->tabu_list[i]);
        }
        free(tabu->tabu_list);
    }

    free(tabu);
}

static int best_2opt_not_tabu(instance *inst, solution *sol, solution *best_sol, tabuList *tabu)
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

            if ((!is_tabu && delta < best_delta && delta != 0) || (is_tabu && new_cost + EPSILON < best_sol->cost))
            {
                if (is_tabu && VERBOSE >= 50)
                    printf("[ASPIRATION] %-34s %10.2f\n", "Accepted move with cost:", new_cost);

                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if (best_i == -1 || best_j == -1)
        return -1;

    reverse_path_segment(best_i, best_j, sol);

    if (evaluate_path_cost(inst, sol))
        return -1;

    add_tabu_move(tabu, best_i, best_j);
    return 0;
}

static tabuList *init_tabu_list(instance *inst)
{
    tabuList *tabu = (tabuList *)malloc(sizeof(tabuList));
    if (!tabu)
    {
        print_error("Failed to allocate tabuList");
        return NULL;
    }

    tabu->min_tenure = inst->tabu_min == 0 ? inst->nnodes / 10 : inst->tabu_min;
    tabu->max_tenure = inst->tabu_max == 0 ? inst->nnodes + inst->nnodes / 10 : inst->tabu_max;
    tabu->tenure = inst->tabu_tenure == 0 ? (tabu->min_tenure + tabu->max_tenure) / 2 : inst->tabu_tenure;

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
    tabu->iterations_without_improvement = 0;

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
    FILE *log_file = NULL;
    if (VERBOSE >= 30)
    {
        log_file = fopen("logs/vns_log.csv", "w");
        if (!log_file)
        {
            print_error("Failed to open log file for VNS");
            return EXIT_FAILURE;
        }
        fprintf(log_file, "iteration,current_cost,best_cost,kicks\n");
    }

    if (!inst || !sol->tour || inst->nnodes <= 0 || !inst->cost_matrix)
        return EXIT_FAILURE;

    double starting_time_vns = second();
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

    if (sol->initialized) // Check if already initialized
    {
        current_sol->cost = sol->cost; // Use the existing cost
        memcpy(current_sol->tour, sol->tour, sizeof(int) * (inst->nnodes + 1));
    }
    else
    {
        // Otherwise, it generates a random path
        if (generate_random_path(inst, current_sol))
        {
            print_error("Random path failed in vns");
            free(current_sol->tour);
            free(current_sol);
            free(best_sol->tour);
            free(best_sol);
            return EXIT_FAILURE;
        }
    }

    memcpy(best_sol->tour, current_sol->tour, sizeof(int) * (inst->nnodes + 1));
    best_sol->cost = current_sol->cost;

    int iter = 0;
    int last_logged_best = -1;
    int no_improvement_count = 0;

    while (!check_time(inst, starting_time_vns))
    {
        if (apply_two_opt(inst, current_sol))
            print_error("2-opt failed in VNS");

        if (current_sol->cost < best_sol->cost)
        {
            best_sol->cost = current_sol->cost;
            memcpy(best_sol->tour, current_sol->tour, sizeof(int) * (inst->nnodes + 1));

            if (VERBOSE >= 10 && iter != last_logged_best)
            {
                printf("[VNS] %-40s%11d\n", "New best found at iter:", iter);
                printf("[VNS] %-40s%11.2f\n", "New best cost:", best_sol->cost);
                last_logged_best = iter;
            }

            no_improvement_count = 0;
        }
        else
        {
            no_improvement_count++;
        }

        int kicks;
        if (inst->vns_jumps > 0)
        {
            kicks = inst->vns_jumps;
        }
        else
        {
            int base_kicks = (rand() % (inst->vns_kmax - inst->vns_kmin + 1)) + inst->vns_kmin;
            kicks = (int)(base_kicks * inst->vns_learning_rate);

            kicks += no_improvement_count / 10;
            if (kicks > inst->nnodes / 2)
                kicks = inst->nnodes / 2;
        }

        if (VERBOSE >= 70)
            printf("[VNS] Applying %d 3-opt kicks for diversification at iter %d\n", kicks, iter);

        for (int i = 0; i < kicks; i++)
        {
            if (apply_three_opt(inst, current_sol))
                print_error("3-opt failed in VNS");
        }

        evaluate_path_cost(inst, current_sol);

        if (log_file)
        {
            fprintf(log_file, "%d,%.6f,%.6f,%d\n", iter, current_sol->cost, best_sol->cost, kicks);
        }

        iter++;
    }

    if (VERBOSE >= 20)
        printf("[VNS] Total elapsed time: %.3fs\n", second() - starting_time_vns);

    if (log_file)
        fclose(log_file);

    memcpy(sol->tour, best_sol->tour, sizeof(int) * (inst->nnodes + 1));
    sol->cost = best_sol->cost;
    sol->initialized = 1; // Mark as initiliazed

    free(current_sol->tour);
    free(current_sol);
    free(best_sol->tour);
    free(best_sol);

    return EXIT_SUCCESS;
}

int apply_heuristic_tabu(instance *inst, solution *sol)
{
    FILE *log_file = NULL;
    if (VERBOSE >= 30)
    {
        log_file = fopen("logs/tabu_log.csv", "w");
        if (!log_file)
        {
            print_error("Failed to open log file for Tabu Search");
            return EXIT_FAILURE;
        }
        fprintf(log_file, "iteration,current_cost,best_cost,tenure\n");
    }

    if (!inst || !sol)
        return EXIT_FAILURE;

    double starting_time = second();

    int nnodes = inst->nnodes;
    tabuList *tabu = init_tabu_list(inst);
    if (!tabu)
        return EXIT_FAILURE;

    if (!sol->initialized)
    {
        if (generate_random_path(inst, sol))
        {
            print_error("Random path failed in tabu");
            free_tabu(tabu);
            return EXIT_FAILURE;
        }
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
    int no_improve_limit = inst->tabu_noimprove == 0 ? NO_IMPROVE_LIMIT : inst->tabu_noimprove;

    while (!check_time(inst, starting_time))
    {
        int did_log = 0;

        if (best_2opt_not_tabu(inst, current_sol, best_sol, tabu) == -1)
            print_error("No valid 2-opt move found");

        if (current_sol->cost < best_sol->cost)
        {
            memcpy(best_sol->tour, current_sol->tour, sizeof(int) * (nnodes + 1));
            best_sol->cost = current_sol->cost;
            tabu->iterations_without_improvement = 0;

            if (VERBOSE >= 10)
            {
                did_log = 1;
                printf("[TABU] %-40s%11d\n", "New best found at iter:", iter);
                printf("[TABU] %-40s%11.2f\n", "New best cost:", best_sol->cost);
            }
        }
        else
        {
            tabu->iterations_without_improvement++;
        }

        // Diversification trigger
        if (tabu->iterations_without_improvement > no_improve_limit)
        {
            if (VERBOSE >= 40)
            {
                did_log = 1;
                printf("[TABU] %-40s%11d\n", "Diversification triggered at:", iter);
            }

            int i = rand() % (nnodes - 2) + 1;
            int j = i + 1 + rand() % (nnodes - i - 1);
            reverse_path_segment(i, j, current_sol);

            if (evaluate_path_cost(inst, current_sol))
                print_error("Failed to evaluate path after diversification");

            // Increase tenure
            int new_tenure = tabu->tenure + 5;
            if (new_tenure > tabu->max_tenure)
                new_tenure = tabu->max_tenure;

            resize_tabu_list(tabu, new_tenure);

            if (VERBOSE >= 60)
            {
                did_log = 1;
                printf("[TABU] %-40s%11d\n", "Increased tenure to:", new_tenure);
            }

            tabu->iterations_without_improvement = 0;
        }
        else if (tabu->iterations_without_improvement == 0 && tabu->tenure > tabu->min_tenure)
        {
            int new_tenure = tabu->tenure - 1;
            if (new_tenure < tabu->min_tenure)
                new_tenure = tabu->min_tenure;

            resize_tabu_list(tabu, new_tenure);

            if (VERBOSE >= 60)
            {
                did_log = 1;
                printf("[TABU] %-40s%11d\n", "Decreased tenure to:", new_tenure);
            }
        }

        if (did_log && VERBOSE >= 10)
            printf("%s\n", DIVIDER);

        iter++;

        if (log_file)
        {
            fprintf(log_file, "%d,%.6f,%.6f,%d\n", iter, current_sol->cost, best_sol->cost, tabu->tenure);
        }
    }

    if (log_file)
        fclose(log_file);

    if (VERBOSE >= 20)
        printf("[TABU] Search terminated due to time limit.\n");

    if (VERBOSE >= 10)
    {
        printf("[TABU] %-40s%11.2f\n", "Final best cost:", best_sol->cost);
        printf("[TABU] %-40s%11.3fs\n", "Total elapsed time:", second() - starting_time);
    }

    memcpy(sol->tour, best_sol->tour, sizeof(int) * (nnodes + 1));
    sol->cost = best_sol->cost;
    sol->initialized = 1;

    free_tabu(tabu);
    free(best_sol->tour);
    free(best_sol);
    free(current_sol->tour);
    free(current_sol);

    return EXIT_SUCCESS;
}
