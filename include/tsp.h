#ifndef TSP_H_
#define TSP_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #include <cplex.h>
#include <pthread.h>

/*
 * printing level:
 *     - 10 only incumbent,
 *     - 20 little output
 *     - 50-60 good,
 *     - 70 verbose,
 *     - 100 cplex log
 */
#define VERBOSE 50

#define EPSILON 1e-9 // very small numerical tolerance
#define XSMALL 1e-5	 // tolerance used to decide ingerality of 0-1 var.s
#define CPX_INFBOUND 1.0E+20

#define MAX_X 10000
#define MAX_Y 10000

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

/**
 * Data structure to define a solution
 */
typedef struct
{
	int *tour;
	double cost;
} solution;

typedef struct
{

	// input data
	int nnodes;
	double *xcoord;
	double *ycoord;

	// computed data
	double starting_time; // starting time
	double **cost_matrix; // cost matrix

	// parameters
	char input_file[1000]; // input file
	char method[20];	   // method to be used
	double timelimit;	   // overall time limit, in sec.s
	int randomseed;

} instance;

void init(instance *inst);
int load_instance(instance *inst);
void print_error(const char *err_message);
void allocate_buffers(instance *tsp);

int parse_command_line(int argc, char **argv, instance *inst);
void free_instance(instance *inst, solution *sol);

void generate_random_nodes(instance *inst, int nnodes, int seed);
int random_path(solution *sol, int nnodes, int seed);
void print_path(const instance *inst, const solution *sol);
void print_nodes(instance *inst);

int check_time(const instance *inst);

int write_path_file(const instance *inst, const solution *sol, const char *filename);

int cost_path(const instance *inst, solution *sol);

int two_opt(const instance *inst, solution *sol);
double delta(int i, int j, const solution *sol, const instance *inst);
void swap_path(int i, int j, solution *sol);

int tsp_compute_costs(instance *tsp);

int run_method(instance *inst, solution *sol);

static void generate_3opt_positions(instance *tsp, int *positions);

int tsp_3opt_solution(instance *tsp, solution *sol);

static void tsp_3opt_swap(instance *tsp, solution *current_sol, int *temp_tour, int *positions);

/*
//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; }
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; }
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; }
*/

#endif /* TSP_H_ */
