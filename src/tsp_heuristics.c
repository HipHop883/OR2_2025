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
        if (generate_random_path(inst, current_sol) != EXIT_SUCCESS)
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

    int min_kicks = inst->vns_kmin;
    int max_kicks = inst->vns_kmax;

    // Gnuplot setup
    FILE *gp = NULL;
    if (inst->plot == 0)
    {
        gp = popen("gnuplot -persist", "w");
        if (gp)
        {
            char filename[256];
            sprintf(filename, "plot/TSP_Heuristic_%s.png", inst->method);

            fprintf(gp, "set terminal png\n");
            fprintf(gp, "set output '%s'\n", filename);
            fprintf(gp, "set title 'TSP VNS Solution Cost vs Iteration'\n");
            fprintf(gp, "set xlabel 'Iteration'\n");
            fprintf(gp, "set ylabel 'Solution Cost'\n");
            fprintf(gp, "plot '-' with lines title 'Costs', '-' with lines lc rgb 'red' title 'Best costs'\n");
        }
    }

    int iter = 0;
    int best_iter = 0;
    double best_cost_local = best_sol->cost;

    // Array to save red line
    int *best_iters = NULL;
    double *best_costs = NULL;
    int best_points_count = 0;
    int best_points_capacity = 0;

    if (inst->plot == 0 && gp)
    {
        best_points_capacity = 100; // Initial Capacity
        best_iters = (int *)malloc(best_points_capacity * sizeof(int));
        best_costs = (double *)malloc(best_points_capacity * sizeof(double));

        if (!best_iters || !best_costs)
        {
            print_error("Memory allocation failed for best_iters/best_costs");
            return EXIT_FAILURE;
        }

        best_iters[best_points_count] = 0;
        best_costs[best_points_count] = best_sol->cost;
        best_points_count++;
    }

    while (!check_time(inst, starting_time_vns))
    {
        if (apply_two_opt(inst, current_sol) != EXIT_SUCCESS)
            print_error("2-opt failed in VNS");

        if (current_sol->cost < best_sol->cost)
        {
            best_sol->cost = current_sol->cost;
            memcpy(best_sol->tour, current_sol->tour, sizeof(int) * (inst->nnodes + 1));
            best_cost_local = best_sol->cost;
            best_iter = iter;
        }

        if (inst->plot == 0 && gp)
        {
            // Check if the it needs to realloc the array
            if (best_points_count >= best_points_capacity)
            {
                best_points_capacity *= 2; // Double the capacity
                int *tmp_iters = (int *)realloc(best_iters, best_points_capacity * sizeof(int));
                double *tmp_costs = (double *)realloc(best_costs, best_points_capacity * sizeof(double));

                if (!tmp_iters || !tmp_costs)
                {
                    print_error("Memory reallocation failed for best_iters/best_costs");
                    free(best_iters);
                    free(best_costs);
                    return EXIT_FAILURE;
                }

                best_iters = tmp_iters;
                best_costs = tmp_costs;
            }

            // Save the new best_cost
            best_iters[best_points_count] = iter;
            best_costs[best_points_count] = best_sol->cost;
            best_points_count++;
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
        }

        for (int i = 0; i < kicks; i++)
        {
            if (apply_three_opt(inst, current_sol) != 0)
                print_error("3-opt failed in VNS");
        }

        evaluate_path_cost(inst, current_sol);

        if (gp)
        {
            fprintf(gp, "%d %lf\n", iter, current_sol->cost);
        }
        else if (inst->plot == 0)
        {
            printf("Gnuplot pipe error\n");
        }
        iter++;
    }

    printf("TIME FOR DO VNS %lf\n", second() - starting_time_vns);

    if (gp)
    {
        fprintf(gp, "e\n");

        // Send all points of the best cost to GnuPlot
        for (int i = 0; i < best_points_count; i++)
        {
            fprintf(gp, "%d %lf\n", best_iters[i], best_costs[i]);
        }
        fprintf(gp, "e\n");
        fflush(gp);
        pclose(gp);
    }
    else if (inst->plot == 0)
    {
        printf("Gnuplot pipe error\n");
    }

    if (VERBOSE >= 20)
        printf("VNS stopped due to time limit.\n");

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
    if (!inst || !sol)
        return EXIT_FAILURE;

    double starting_time = second();

    int nnodes = inst->nnodes;
    tabuList *tabu = init_tabu_list(nnodes);
    if (!tabu)
        return EXIT_FAILURE;

    if (!sol->initialized) // Check if solution is not initialized
    {
        if (generate_random_path(inst, sol) != 0)
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

    while (!check_time(inst, starting_time))
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
    sol->initialized = 1; // Mark as initialized

    // free_tabu(tabu);
    free(best_sol->tour);
    free(best_sol);
    free(current_sol->tour);
    free(current_sol);

    return EXIT_SUCCESS;
}
