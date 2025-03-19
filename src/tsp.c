#include "../include/tsp.h"
#include "chrono.h"

/** 
 * Read input data
 * @param inst instance
 * @return 0 if the input data is read successfully, 1 otherwise
*/
int read_input(instance *inst) {
    if (strcmp(inst->input_file, "NULL") == 0) {
        // No input file
		if (inst->nnodes <= 0) {
            printf("Please provide the number of nodes: ");
            scanf("%d", &inst->nnodes);
        }
        if (inst->randomseed <= 0) {
            printf("Please provide the seed for the random generator: ");
            scanf("%d", &inst->randomseed);
        }
        generate_random_nodes(inst, inst->nnodes, inst->randomseed);

		if (strcmp(inst->method, "NULL") == 0) {
			printf("Please provide the method to be used: ");
			scanf("%s", inst->method);	
		}
        return EXIT_SUCCESS;
    }
    
    FILE *fin = fopen(inst->input_file, "r");
    if (fin == NULL) 
	{
		print_error("File not found");
		return EXIT_FAILURE;
	}

    inst->nnodes = -1;

    char line[100];
    char *par_name;   
	char *token1;
	char *token2;
	
	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION 
	
	//int do_print = ( VERBOSE >= 1000 );

	while ( fgets(line, sizeof(line), fin) != NULL ) 
	{
		if ( VERBOSE >= 2000 ) { printf("%s",line); fflush(NULL); }
		if ( strlen(line) <= 1 ) continue; // skip empty lines
	    par_name = strtok(line, " :");
		if ( VERBOSE >= 3000 ) { printf("parameter \"%s\" ",par_name); fflush(NULL); }

		if ( strncmp(par_name, "NAME", 4) == 0 ) 
		{
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "COMMENT", 7) == 0 ) 
		{
			active_section = 0;   
			token1 = strtok(NULL, "");  
			// if ( VERBOSE >= 10 ) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
			continue;
		}   
		
		if ( strncmp(par_name, "TYPE", 4) == 0 ) 
		{
			token1 = strtok(NULL, " :");  
			if ( strncmp(token1, "TSP",3) != 0 ) 
			{
				print_error(" format error:  only TYPE == TSP implemented so far!!!!!!"); 
				return EXIT_FAILURE;
			}
			active_section = 0;
			continue;
		}
		

		if ( strncmp(par_name, "DIMENSION", 9) == 0 ) 
		{
			if ( inst->nnodes >= 0 ) {
				print_error(" repeated DIMENSION section in input file");
				return EXIT_FAILURE;
			}
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			//if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes); 
			//inst->demand = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double)); 
			inst->cost_matrix = (double **) calloc(inst->nnodes, sizeof(double));
			for (int i = 0; i < inst->nnodes; i++) {
				inst->cost_matrix[i] = (double *) calloc(inst->nnodes, sizeof(double));
			}

			active_section = 0;  
			continue;
		}

		if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 ) 
		{
			token1 = strtok(NULL, " :");
			if ( strncmp(token1, "ATT", 3) != 0 ) 
			{
				print_error(" format error:  only EDGE_WEIGHT_TYPE == ATT implemented so far!!!!!!"); 
				return EXIT_FAILURE;
			}
			active_section = 0;
			continue;
		}            
		
		if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 ) 
		{
			if ( inst->nnodes <= 0 ) 
			{
				print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
				return EXIT_FAILURE;
			}
			active_section = 1;   
			continue;
		}
		/*
		if ( strncmp(par_name, "DEMAND_SECTION", 14) == 0 ) 
		{
			if ( inst->nnodes <= 0 ) {
			print_error(" ... DIMENSION section should appear before DEMAND_SECTION section");
			return EXIT_FAILURE;
			}
			active_section = 2;
			continue;
		}  
        
		if ( strncmp(par_name, "DEPOT_SECTION", 13) == 0 )  
		{
			if ( inst->depot >= 0 ) {
			print_error(" ... DEPOT_SECTION repeated??");
			return EXIT_FAILURE;
			}
			active_section = 3;   
			continue;
		}
		*/

		if ( strncmp(par_name, "EOF", 3) == 0 ) 
		{
			active_section = 0;
			break;
		}
		
			
		if ( active_section == 1 ) // within NODE_COORD_SECTION
		{
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) 
			{
				print_error(" ... unknown node in NODE_COORD_SECTION section");     
				return EXIT_FAILURE;
			}
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			//if ( do_print ) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i+1, inst->xcoord[i], inst->ycoord[i]); 
			continue;
		}    
		/*
		if ( active_section == 2 ) // within DEMAND_SECTION
		{
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) {
			print_error(" ... unknown node in NODE_COORD_SECTION section");     
			return EXIT_FAILURE;
			}	
			token1 = strtok(NULL, " :,");
			inst->demand[i] = atof(token1);
			if ( do_print ) printf(" ... node %4d has demand %10.5lf\n", i+1, inst->demand[i]); 
			continue;
		}  

		if ( active_section == 3 ) // within DEPOT_SECTION
		{
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) continue;
			if ( inst->depot >= 0 ) 
			{
			print_error(" ... multiple depots not supported in DEPOT_SECTION");     
			return EXIT_FAILURE;
			}	
			inst->depot = i;
			if ( do_print ) printf(" ... depot node %d\n", inst->depot+1); 
			continue;
		}  
		*/
		printf(" final active section %d\n", active_section);
		print_error(" ... wrong format for the current simplified parser!!!!!!!!!");     
		return EXIT_FAILURE;
		    
	}                

	// Compute cost matrix
	if (tsp_compute_costs(inst) != 0) {
		print_error("Error computing cost matrix");
		return EXIT_FAILURE;
	}

	fclose(fin);
	return EXIT_SUCCESS;
}

/**
 * Print error message
 * @param err_message error message
 * @return void
*/
void print_error(const char *err_message) { fprintf(stderr, "\n\n ERROR: %s \n\n", err_message);}

/** 
 * Parse command line
 * @param argc number of arguments
 * @param argv arguments
 * @param inst instance
 * @return 0 if the command line is parsed successfully, 1 otherwise
*/
int parse_command_line(int argc, char** argv, instance *inst) 
{ 
	
	if ( VERBOSE >= 100 ) printf(" running %s with %d parameters \n", argv[0], argc-1); 
	
	// default   
	strcpy(inst->input_file, "NULL");
	inst->randomseed = 0;     
    inst->timelimit = CPX_INFBOUND;
	inst->nnodes = 0;
	strcpy(inst->method, "NULL");
	inst->cost_matrix = NULL;

    int help = 0;
    if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{ 
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
        if ( strcmp(argv[i],"-nnodes") == 0 && strcmp(inst->input_file, "NULL")==0) 
			{ inst->nnodes = abs(atoi(argv[++i])); continue; } 											// number of nodes
		if ( strcmp(argv[i],"-method") == 0 ) { strcpy(inst->method, argv[++i]); continue; } 			// method
		if ( strcmp(argv[i],"-m") == 0 ) { strcpy(inst->method, argv[++i]); continue; } 				// method
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
    }

	if ( help || (VERBOSE >= 10) )		// print current parameters
	{
		printf("\n\navailable parameters (vers. 03-march-2025) --------------------------------------------------\n");
		printf("-file %s\n", inst->input_file); 
		printf("-time_limit %lf\n", inst->timelimit); 
		printf("-seed %d\n", inst->randomseed); 
		printf("-method %s\n", inst->method);
		printf("\nenter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}  
	
	if ( help ) exit(1);
/*
	else if ( strcmp(inst->input_file, "NULL") == 0) { // Without an input file
        if (inst->nnodes <= 0) {
            printf("Please provide the number of nodes: ");
            scanf("%d", &inst->nnodes);
        }
        if (inst->randomseed <= 0) {
            printf("Please provide the seed for the random generator: ");
            scanf("%d", &inst->randomseed);
        }
        generate_random_nodes(inst, inst->nnodes, inst->randomseed);

		if (strcmp(inst->method, "NULL") == 0) {
			printf("Please provide the method to be used: ");
			scanf("%s", inst->method);	
		}
    }
*/   
	return EXIT_SUCCESS;
} 

/**
 * Generate random nodes
 * @param inst instance
 * @param nnodes number of nodes
 * @param seed seed
 * @return void
 */
void generate_random_nodes(instance *inst, int nnodes, int seed) {
	if (VERBOSE >=50) printf("Generating random nodes...\n");

    inst->nnodes = nnodes;
    inst->xcoord = (double *) calloc(nnodes, sizeof(double));
    inst->ycoord = (double *) calloc(nnodes, sizeof(double));
    
    srand(seed);
    for (int i = 0; i < nnodes; i++) {
        inst->xcoord[i] = ((double) rand() / RAND_MAX) * MAX_X;
        inst->ycoord[i] = ((double) rand() / RAND_MAX) * MAX_Y;
    }

	if (VERBOSE >=50) printf("Random nodes generated\n");
}

/**
 * Generate random nodes
 * @param inst instance
 * @return void
 
void generate_random_nodes(instance *inst) {
	srand(7654321 + inst->randomseed);
	
	inst->nnodes = (int) (rand() % RAND_MAX)*1000;
	inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double));
	inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));
	for (int i = 0; i < inst->nnodes; i++) {
		inst->xcoord[i] = ((double) rand() / RAND_MAX) * MAX_X;
		inst->ycoord[i] = ((double) rand() / RAND_MAX) * MAX_Y;
	}
}
*/

/** 
 * Free instance
 * @param inst instance
 * @return void
*/
void free_instance(instance *inst, solution *sol)
{     
	free(inst->xcoord);
	inst->xcoord = NULL;

	free(inst->ycoord);
	inst->ycoord = NULL;

	for (int i = 0; i < inst->nnodes; i++) {
        free(inst->cost_matrix[i]);
	}
	free(inst->cost_matrix);
	inst->cost_matrix = NULL;

	free(sol->tour);
	sol->tour = NULL;
}

/**
 * Generate a random path using the Fisher-Yates shuffle algorithm
 * @param sol path of the solution
 * @param nnodes number of nodes
 * @param seed seed
 * @return 0 if the random path is generated successfully, 1 otherwise
 */
int random_path(solution *sol, int nnodes, int seed) {
    int nodes[nnodes];
    for (int i = 0; i < nnodes; i++) {
        nodes[i] = i;
    }

    srand(seed);
    for (int i = nnodes - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        // Swap nodes[i] with nodes[j]
        int temp = nodes[i];
        nodes[i] = nodes[j];
        nodes[j] = temp;
    }

    for (int i = 0; i < nnodes; i++) {
        sol->tour[i] = nodes[i];
    }
    sol->tour[nnodes] = sol->tour[0]; // Comes back to the initial node
	return EXIT_SUCCESS;
}

/** 
 * Print path
 * @param inst instance
 * @param sol solution path
 * @return void
*/
void print_path(const instance *inst, const solution *sol) {
    for (int i = 0; i <= inst->nnodes; i++) 
	{
        int node = sol->tour[i];
        printf("(%d, %d) ", (int) inst->xcoord[node], (int)inst->ycoord[node]);
    }
    printf("\n");
}

/**
 * Print nodes
 * @param inst instance
 * @return void
 */
void print_nodes(instance *inst) {
	for (int i = 0; i < inst->nnodes; i++) {
		printf("Node %d: (%d, %d)\n", i, (int)inst->xcoord[i], (int)inst->ycoord[i]);
	}
}

/**
 * Check time
 * @param inst instance
 * @param start_time start time
 * @return 0 if the time limit is not reached, 1 otherwise
 */
int check_time(const instance *inst) {
	double current_time = second();
	if (current_time - inst->starting_time > inst->timelimit) {
		print_error("Time limit reached");
		exit(1);
	}
	return EXIT_SUCCESS;
}

/**
 * Write path file
 * @param inst instance
 * @param sol solution path
 * @param filename filename
 * @return 0 if the path file is written successfully, 1 otherwise
 */
int write_path_file(const instance *inst, const solution *sol, const char *filename) {
    FILE *fout = fopen(filename, "w");
    if (fout == NULL) 
	{
		print_error("File not found");
		return EXIT_FAILURE;
	}

    for (int i = 0; i <= inst->nnodes; i++) {
        int node = sol->tour[i];
        fprintf(fout, "%d %d\n", (int)inst->xcoord[node], (int)inst->ycoord[node]);
    }

    fclose(fout);
	return EXIT_SUCCESS;
}

/**
 * cost of the edge between nodes i and j
 * @param i node i
 * @param j node j
 * @param inst instance
 * @return cost of the edge between nodes i and j
 */

double cost(int i, int j, instance *inst) {																		// I GUESS TO BE REMOVED
	return sqrt(pow(inst->xcoord[i] - inst->xcoord[j], 2) + pow(inst->ycoord[i] - inst->ycoord[j], 2));
}

/**
 * cost of the path
 * @param inst instance
 * @param sol solution path
 * @return 0 if the cost of the path is calculated successfully, 1 otherwise
 */

int cost_path(const instance *inst, solution *sol) {
    double cost_p = 0;
    for (int i = 0; i < inst->nnodes - 1; i++) {
        //double edge_cost = cost(inst->best_sol[i], inst->best_sol[i + 1], inst);															// TO BE REMOVED
        //double edge_cost = inst->cost_matrix[flatten_coords(inst->best_sol[i], inst->best_sol[i + 1], inst->nnodes)]; 					// TO BE REMOVED
        double edge_cost = inst->cost_matrix[sol->tour[i]][sol->tour[i + 1]];
        
		if (VERBOSE >= 70) printf("Cost from node %d to node %d: %lf\n", sol->tour[i], sol->tour[i + 1], edge_cost);
        cost_p += edge_cost;
    }

    // Add the cost to return to the starting node
    //double return_cost = cost(inst->best_sol[inst->nnodes - 1], inst->best_sol[0], inst); 											// TO BE REMOVED
    //double return_cost = inst->cost_matrix[flatten_coords(inst->best_sol[inst->nnodes - 1], inst->best_sol[0], inst->nnodes)]; 		// TO BE REMOVED
    double return_cost = inst->cost_matrix[sol->tour[inst->nnodes - 1]][sol->tour[0]];

    if (VERBOSE >= 70) printf("Cost from node %d to node %d: %lf\n", sol->tour[inst->nnodes - 1], sol->tour[0], return_cost);
    cost_p += return_cost;
    printf("Total cost: %lf\n", cost_p);
	sol->cost = cost_p;			// update the cost of the solution

	return EXIT_SUCCESS;
}

/**
 * Nearest neighbor heuristic
 * @param inst instance
 * @param sol solution path
 * @return 0 if the nearest neighbor heuristic is applied successfully, 1 otherwise
*/ 
int nearest_neighbor(const instance *inst, solution *sol) {
    double time = second();
    if (VERBOSE >= 50) {
        printf("Applying nearest neighbor heuristic...\n");
        printf("Starting time: %lf\n", time);
    }

    int nnodes = inst->nnodes;                            // number of nodes
    int *visited = (int *)calloc(nnodes, sizeof(int));    // visited nodes
    int current = 0;                                      // current node
    visited[current] = 1;                                 // mark the current node as visited
    sol->tour[0] = current;                          	  // start from the current node

    for (int i = 1; i < nnodes; i++) {
        double min_cost = CPX_INFBOUND;
        int next = -1;
        for (int j = 0; j < nnodes; j++) {  // find the nearest neighbor
            if (visited[j] == 0) {
                //double c = inst->cost_matrix[flatten_coords(current, j, nnodes)];
                double c = inst->cost_matrix[current][j];
                if (c < min_cost) {  // update the nearest neighbor
                    min_cost = c;    // update the minimum cost
                    next = j;        // update the next node
                }
            }
        }
        visited[next] = 1;           // mark the next node as visited
        sol->tour[i] = next;    // update the best solution
        current = next;              // update the current node
    }
    sol->tour[nnodes] = sol->tour[0];  // return to the initial node
    free(visited);

    if (VERBOSE >= 50) {
        printf("Nearest neighbor heuristic applied\n");
        printf("Elapsed time: %lf\n", second() - time);
    }

	return EXIT_SUCCESS;
}

/**
 * Two-opt heuristic
 * @param inst instance
 * @param sol solution path
 * @return 0 if the two-opt heuristic is applied successfully, 1 otherwise
 */
int two_opt(const instance *inst, solution *sol) {
    double time = second();
    if (VERBOSE >= 50) {
        printf("Applying two-opt heuristic...\n");
        printf("Starting time: %lf\n", time);
    }
    int nnodes = inst->nnodes;
    int improved = 1;  		// flag to indicate if the path has been improved
    while (improved) {
        improved = 0;  		// reset the flag
        for (int i = 0; i < nnodes - 1; i++) {
            for (int j = i + 1; j < nnodes; j++) {
                if (delta(i, j, sol, inst) < 0) {  // if the path is improved
                    swap_path(i, j, sol);          // swap nodes i and j
                    improved = 1;                  // set the flag
                }
            }
        }
    }

    if (VERBOSE >= 50) {
        printf("Two-opt heuristic applied\n");
        printf("Elapsed time: %lf\n", second() - time);
    }
	
	return EXIT_SUCCESS;
}

/**
 * Calculate the change in the cost of the path when swapping nodes i and j
 * @param i node i
 * @param j node j
 * @param sol solution path
 * @param inst instance
 * @return change in the cost of the path when swapping nodes i and j
 */

double delta(int i, int j, const solution *sol, const instance *inst) {
	/*return (inst->cost_matrix[flatten_coords(path[i+1], path[j+1], inst->nnodes)] + inst->cost_matrix[flatten_coords(path[i], path[j], inst->nnodes)] 				// TO BE REMOVED
		-(inst->cost_matrix[flatten_coords(path[i], path[i+1], inst->nnodes)] + inst->cost_matrix[flatten_coords(path[j], path[j+1], inst->nnodes)]));
	*/
	return (inst->cost_matrix[sol->tour[i+1]][sol->tour[j+1]] + inst->cost_matrix[sol->tour[i]][sol->tour[j]] 
		-(inst->cost_matrix[sol->tour[i]][sol->tour[i+1]] + inst->cost_matrix[sol->tour[j]][sol->tour[j+1]]));
}

/**	
 * Swap nodes i and j in the path
 * @param i node i
 * @param j node j
 * @param sol solution path
 * @return void
 */
void swap_path(int i, int j, solution *sol) {
	while (++i < j) {
		double temp = sol->tour[i];
		sol->tour[i] = sol->tour[j];
		sol->tour[j] = temp;
		j--;
	}
}



/**
 * Compute the cost matrix
 * @param tsp instance
 * @return 0 if the cost matrix is computed successfully, 1 otherwise
 */
int tsp_compute_costs(instance *tsp)
{
    if (tsp->cost_matrix == NULL || tsp->xcoord == NULL || tsp->ycoord == NULL)
	{
		return -1;
	}

    for (int i = 0; i < tsp->nnodes; i++)
    {
        for (int j = 0; j < tsp->nnodes; j++)
        {
            double deltax = tsp->xcoord[i] - tsp->xcoord[j];
            double deltay = tsp->ycoord[i] - tsp->ycoord[j];
            double dist = sqrt(deltax * deltax + deltay * deltay);

            //tsp->cost_matrix[flatten_coords(i, j, tsp->nnodes)] = dist;						// TO BE REMOVED
            tsp->cost_matrix[i][j] = dist;
        }
    }

    return 0;
}

/**
 * Run the method
 * @param inst instance
 * @param sol solution path
 * @return 0 if the method is run successfully, 1 otherwise
 */
int run_method(instance *inst, solution *sol)
{
	if(strcmp(inst->method, "n_n") == 0) {
		if(nearest_neighbor(inst, sol)) 			// Nearest neighbor heuristic
		{
			print_error("Error applying nearest neighbor heuristic");
			return EXIT_FAILURE;
		}
		if(check_time(inst)) 
		{
			return EXIT_FAILURE;
		}
	} else if(strcmp(inst->method, "n_n+two_opt") == 0) {
		if(nearest_neighbor(inst, sol))			// Nearest neighbor heuristic with two_opt
		{
			print_error("Error applying nearest neighbor heuristic");
			return EXIT_FAILURE;
		}
		if(two_opt(inst, sol)) 
		{
			print_error("Error applying 2-opt");
			return EXIT_FAILURE;
		}
		if (check_time(inst))
		{
			return EXIT_FAILURE;
		}
	} else if(strcmp(inst->method, "random+two_opt") == 0) {			// Random path with two_opt
		if(random_path(sol, inst->nnodes, inst->randomseed))
		{
			print_error("Error generating random path");
			return EXIT_FAILURE;
		}	
		if(two_opt(inst, sol))
		{
			print_error("Error applying 2-opt");
			return EXIT_FAILURE;
		}
		if (check_time(inst))
		{
			return EXIT_FAILURE;
		}
	} else if(strcmp(inst->method, "random") == 0) {
		random_path(sol, inst->nnodes, inst->randomseed);	// Random path
		if (check_time(inst)) 
		{
			return EXIT_FAILURE;
		}
	} else {
		print_error("USAGE: -method [n_n|n_n+two_opt|random+two_opt|random]");
		return EXIT_FAILURE;
	}

	// Print the cost of the best solution
	cost_path(inst, sol);
	return EXIT_SUCCESS;
}


