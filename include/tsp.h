#ifndef TSP_H_  

#define TSP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 

//#include <cplex.h>  
#include <pthread.h>  

#define VERBOSE				    50		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define CPX_INFBOUND  		  1.0E+20
//#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ
                                 
#define MAX_X 100
#define MAX_Y 1000

#define flatten_coords(x, y, N) (x * N + y)

//data structures  

typedef struct {   
	
	//input data
	int nnodes; 	   
	double *xcoord;
	double *ycoord;

	// parameters 
 
	int randomseed;
	double timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file
	char method[20];						// method to be used
	double starting_time;					// starting time
	double **cost_matrix;					// cost matrix
	
	//global variables
	int *best_sol;						// best sol. available
	
} instance; 

/**
 * Data structure to define a solution
*/

typedef struct
{
    int *tour;
    double cost;

} solution;

void read_input(instance *inst);
void print_error(const char *err_message);

void parse_command_line(int argc, char** argv, instance *inst);
void free_instance(instance *inst);

void generate_random_nodes(instance *inst, int nnodes, int seed);
void random_path(int *path, int nnodes, int seed);
void print_path(instance *inst, int *path, int nnodes);
void print_nodes(instance *inst);

void check_time(instance *inst);

void write_path_file(instance *inst, int *best_sol, const char *filename);

double cost(int i, int j, instance *inst);
double cost_path(instance *inst);
void nearest_neighbor(instance *inst);

void two_opt(instance *inst);
double delta(int i, int j, int *path, instance *inst);
void swap_path(int i, int j, int *path);

int tsp_compute_costs(instance *tsp);

/*
//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; } 
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; } 
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; } 
*/

#endif   /* TSP_H_ */ 
