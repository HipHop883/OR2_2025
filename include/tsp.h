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
#define VERBOSE 70

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
	int *tour;
	double cost;
	int initialized;
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
	char method[40];
	double timelimit;
	int randomseed;
	int plot;	  // to plot the solution
	int nmethods; // number of methods
	char csv_filename[256];

	// Internal state
	double **cost_matrix;
	double starting_time;
	int edge_weight; // 0 for ATT, 1 for EUC_2D

	// Greedy params
	int greedy_starts;

	// VNS params
	int vns_kmin;
	int vns_kmax;
	double vns_learning_rate;
	int vns_jumps;

	// TABU params
	int tabu_tenure;
	int tabu_min;
	int tabu_max;
	int tabu_noimprove;

	// CPLEX model metadata
	int ncols;
	
	// Hard Fixing percentage
	double hard_fixing_percentage;
} instance;

/**
 * Types of 3-opt reconnection moves.
 * Defines how the three segments (A, B, C) are reversed or not.
 */
enum ThreeOptType
{
	TYPE_0, ///< Identity
	TYPE_1, ///< Reverse A
	TYPE_2, ///< Reverse B
	TYPE_3, ///< Reverse C
	TYPE_4, ///< Reverse A and B
	TYPE_5, ///< Reverse A and C
	TYPE_6	///< Reverse B and C
};

void init(instance *inst);
int parse_command_line(int argc, char **argv, instance *inst);

int load_instance(instance *inst);
void allocate_buffers(instance *inst);
void free_instance(instance *inst);
void free_sol(solution *sol);

void generate_random_nodes(instance *inst, int nnodes);
int generate_random_path(const instance *inst, solution *sol);
int check_time(const instance *inst, double starting_time);

void print_solution_path(const instance *inst, const solution *sol);
void print_node_coordinates(instance *inst);
int write_path_to_file(const instance *inst, const solution *sol, const char *filename);

int evaluate_path_cost(const instance *inst, solution *sol);

int apply_two_opt(const instance *inst, solution *sol);
double path_cost_delta(int i, int j, const solution *sol, const instance *inst);
void reverse_path_segment(int i, int j, solution *sol);

int tsp_compute_costs(instance *tsp);
int execute_selected_method(instance *inst, solution *sol);
double get_method_weight(const char *method_name);


int apply_three_opt(instance *tsp, solution *sol);

void print_error(const char *err_message);

#endif /* TSP_H_ */
