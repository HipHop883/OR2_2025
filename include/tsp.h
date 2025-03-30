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
 * Structure representing a TSP solution (path + cost).
 */
typedef struct
{
	int *tour;	 // sequence of nodes
	double cost; // total path cost
} solution;

/**
 * Structure representing a TSP problem instance.
 */
typedef struct
{
	// Input
	int nnodes;
	double *xcoord;
	double *ycoord;

	// Parameters
	char input_file[1000];
	char method[20];
	double timelimit;
	int randomseed;

	// Internal state
	double **cost_matrix;
	double starting_time;
} instance;

void init(instance *inst);
int parse_command_line(int argc, char **argv, instance *inst);

int load_instance(instance *inst);
void allocate_buffers(instance *inst);
void free_instance(instance *inst, solution *sol);

void generate_random_nodes(instance *inst, int nnodes, int seed);
int generate_random_path(solution *sol, int nnodes, int seed);

int apply_two_opt(const instance *inst, solution *sol);
int apply_three_opt(instance *inst, solution *sol); // wrapper for 3-opt move

double path_cost_delta(int i, int j, const solution *sol, const instance *inst);
void reverse_path_segment(int i, int j, solution *sol);

int tsp_compute_costs(instance *inst);
int evaluate_path_cost(const instance *inst, solution *sol);
int execute_selected_method(instance *inst, solution *sol);

int check_time(const instance *inst);

void print_solution_path(const instance *inst, const solution *sol);
void print_node_coordinates(instance *inst);
int write_path_to_file(const instance *inst, const solution *sol, const char *filename);

void print_error(const char *err_message);

#endif /* TSP_H_ */
