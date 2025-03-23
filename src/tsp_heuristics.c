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

    tabuList *tabu = (tabuList *)malloc(sizeof(tabuList));
    if (!tabu) {
        fprintf(stderr, "memory allocation failed for tabu struct.\n");
        return 1;
    }
    tabu->tenure = MIN_TENURE;                                            //tenure initialized to the minimum value
    tabu->tabu_list = (int **) calloc(tabu->tenure, sizeof(int));         //tabu list allocated
    for(int i = 0; i < tabu->tenure; i++)
    {
        tabu->tabu_list[i] = (int *) calloc(2, sizeof(int));
        if (tabu->tabu_list[i] == NULL)
        {
            fprintf(stderr, "memory allocation failed for tabu list.\n");
            // Free allocated memory
            for (int j = 0; j < i; j++) {
                free(tabu->tabu_list[j]);
            }
            free(tabu->tabu_list);
            free(tabu);
            return 1;
        }
        // Initialize the tabu list
        tabu->tabu_list[i][0] = -1;
        tabu->tabu_list[i][1] = -1;
    }

    // Number of elements in the list
    tabu->size = 0;

    // Get an initial solution using the greedy algorithm.
    if (nearest_neighbor(inst, sol, 0) != 0)
    {
        print_error("Nearest neighbor heuristic failed in tabu");
        free_tabu(tabu);
        return 1;
    }

    if (two_opt(inst, sol) != 0)
    {
        print_error("Two Opt failed in tabu");
        free_tabu(tabu);
        return 1;
    }

    // Copy the initial solution to the best solution
    solution *best_sol = (solution *) malloc(sizeof(solution));
    best_sol->tour = (int *) malloc((nnodes + 1) * sizeof(int));
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
    solution *current_sol = (solution *) malloc(sizeof(solution));
    current_sol->tour = (int *) malloc((nnodes + 1) * sizeof(int));
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
    printf("Initial current cost: %lf\n", current_sol->cost);


    // Initialize the number of iterations
    int iter = 0;    

    while(second() - inst->starting_time < inst->timelimit)
    {
        if (two_opt_random(inst, current_sol, tabu, iter) == -1)
            print_error("Error in two opt random");
        if (cost_path(inst, current_sol))
        {
            print_error("computing cost_path");
            return -1;
        }

        if (current_sol->cost < best_sol->cost)
        {
            memcpy(best_sol->tour, current_sol->tour, (inst->nnodes + 1) * sizeof(int));
            best_sol->cost = current_sol->cost;
        }
        iter++;
    }

    // Copy the best found solution to the output parameters.
    memcpy(sol->tour, best_sol->tour, sizeof(int) * (inst->nnodes + 1));
    sol->cost = best_sol->cost;

    //free memory
    free_tabu(tabu);
    free(best_sol->tour);
    free(best_sol);
    free(current_sol->tour);
    free(current_sol);

    return 0;
}

/** 
 * Apply a 2-opt random move on the current solution.
 * @param tsp TSP instance
 * @param sol solution
 * @param tabu tabu list struct
 * @return 0 if the 2-opt move is applied successfully, 1 if the move is in the tabu list, -1 otherwise
 */
static int two_opt_random(instance *tsp, solution *sol, tabuList *tabu, int iter)
{
    if (!tsp->cost_matrix || tsp->nnodes <= 0)
    {
        print_error("Empty instance");
        return -1;
    }

    if (!tabu->tabu_list)
    {
        print_error("Empty tabu list");
        return -1;
    }

    int i, j;
    int max_attempts = tsp->nnodes * tsp->nnodes;
    int attempts = 0;

    while (attempts < max_attempts) 
    {
        i = rand() % tsp->nnodes;
        j = rand() % tsp->nnodes;

        if (i == j)
        {
            j = (j + 1) % tsp->nnodes;
        }

        if (i > j)
        {
            int temp = i;
            i = j;
            j = temp;
        }
        
        if (!check_tabu_list(tabu, i, j, iter))         //check if the move (i, j) is in the tabu list
        {
            if (tabu->size < tabu->tenure) 
            {
                swap_path(i, j, sol);                   //swap nodes i and j
                tabu->tabu_list[tabu->size][0] = i;     //add the move (i, j) to the tabu list
                tabu->tabu_list[tabu->size][1] = j;
                tabu->size++;
            }

            if (cost_path(tsp, sol))
            {
                print_error("Error computing the (two opt random) cost of the solution");
                return -1;
            }

            return 0;
        }   
        attempts++;
    }

    return -1;     // Can't find a tabu move
}

/**
 * Check if the move (i, j) is in the tabu list.
 * @param tabu_list tabu list
 * @param i node i
 * @param j node j
 * @return 1 if the move (i, j) is in the tabu list, 0 otherwise
 */
int check_tabu_list(tabuList *tabu, int i, int j, int iter) 
{
    for (int k = 0; k < tabu->size; k++) {
        if (tabu->tabu_list[k][0] == i && tabu->tabu_list[k][1] == j) {         // (i, j) is in the tabu list
            return 1;
        }
    }
    if (iter >= MIN_TENURE) 
    {
        int cycle = (iter - MIN_TENURE) % ((MAX_TENURE - MIN_TENURE) * 2); // Cycle of 2*(MAX_TENURE-MIN_TENURE) iterations (half increase + half decrease)

        if (cycle < MAX_TENURE - MIN_TENURE) { // Increase for the first half iterations
            if (tabu->tenure < MAX_TENURE) {
                if (increase_list(tabu) != 0) {
                    print_error("Error increasing the tabu list");
                    return -1; // It is a tabu move
                }
            }
        } else { // Decrease for other MAX_TENURE-MIN_TENURE iterations
            if (tabu->tenure > MIN_TENURE) {
                if (decrease_list(tabu) != 0) {
                    print_error("Error decreasing the tabu list");
                    return -1; // In case of error it is a tabu move
                }
            }
        }
    }

    return 0;
}

/**
 * Increase the tabu list.
 * @param tabu tabu list
 * @return 0 if the tabu list is increased successfully, -1 for errors
 */
int increase_list(tabuList *tabu)
{
    if (tabu->tenure <= 0 || tabu->tabu_list == NULL) {
        print_error("Invalid tabu list");
        return -1;
    }
    tabu->tenure++;
    int **temp;    
    
    //reallocating memory for the tabu list
    temp = (int **) realloc(tabu->tabu_list, tabu->tenure * sizeof(int *));
    
    //check if memory allocation failed
    if(!temp)
    {
        print_error("Memory allocation failed");
        return -1;
    }

    //initialize the new line 
    temp[tabu->tenure - 1] = (int *) calloc(2, sizeof(int));
    if(!temp[tabu->tenure - 1])
    {
        print_error("Memory allocation failed");
        return -1;
    }
    temp[tabu->tenure-1][0] = -1;
    temp[tabu->tenure-1][1] = -1;

    tabu->tabu_list = temp;

    return 0;
}

/**
 * Decrease the tabu list
 * @param tabu tabu list
 * @return 0 if the tabu list is decreased successfully, -1 for errors
 */
int decrease_list(tabuList *tabu)
{
    if (tabu->tenure < MIN_TENURE || tabu->tabu_list == NULL) {
        print_error("Invalid tabu list");
        return -1;
    }

    free(tabu->tabu_list[0]);

    for(int i = 1; i < tabu->tenure; i++)
    {
        tabu->tabu_list[i - 1] = tabu->tabu_list[i];
    }

    tabu->size--;
    tabu->tenure--;

    // Check if tenure == 0
    if (tabu->tenure == 0) {
        free(tabu->tabu_list);
        tabu->tabu_list = NULL;
        print_error("TENURE == 0");
        return -1;
    }

    //reallocating memory for the tabu list
    int **temp = (int **) realloc(tabu->tabu_list, tabu->tenure * sizeof(int *));
    
    //check if memory allocation failed
    if(!temp && tabu->tenure > 0)
    {
        print_error("Memory allocation failed");
        return -1;
    }

    tabu->tabu_list = temp;

    return 0;
}

/**
 * Method to free the tabu list struct
 * @param tabu tabu list struct
 * @return void
 */
void free_tabu(tabuList *tabu) {
    for(int i = 0; i < tabu->tenure; i++) free(tabu->tabu_list[i]);
        free(tabu->tabu_list);
        free(tabu);
}
