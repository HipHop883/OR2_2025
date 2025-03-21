#include <tabu.h>

#define TENURE 10

/**
 * 
 * @param inst instance
 * @param sol solution
 * @return 0 if the solution is found, 1 otherwise
 */
int tabu_method(instance *inst, solution *sol)
{
    if (!inst || !sol)
    {
        return 1;
    }

    int nnodes = inst->nnodes;

    int **tabu_list = (int **) calloc(TENURE, sizeof(int));         //tabu list allocated
    for(int i = 0; i < TENURE; i++)
    {
        tabu_list[i] = (int *) calloc(2, sizeof(int));
    }

    if (tabu_list == NULL)
    {
        fprintf(stderr, "memory allocation failed for tabu list.\n");
        return 1;
    }



    //free memory
    for(int i = 0; i < TENURE; i++) free(tabu_list[i]);
    free(tabu_list);
}

/** 
 * Apply a 2-opt random move on the current solution.
 * @param tsp TSP instance
 * @param sol solution
 * @return 0 if the 2-opt move is applied successfully, 1 otherwise
 */
static int two_opt_random(instance *tsp, solution *sol, int **tabu_list)
{
    if (!tsp->cost_matrix || tsp->nnodes <= 0)
    {
        print_error("Cost matrix not computed or empty instance");
        return 1;
    }

    int i = rand() % tsp->nnodes;
    int j = rand() % tsp->nnodes;

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

    swap_path(i, j, sol);
    if (cost_path(tsp, sol))
    {
        print_error("Error computing the (two opt random) cost of the solution");
        return 1;
    }

    return 0;
}
