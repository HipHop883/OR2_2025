#include "tsp.h"
#include "utils.h"
#include "tsp_heuristics.h"
#include "tsp_greedy.h"
#include "tsp_cplex.h"

/**
 * Load instance data from file or generate randomly
 *
 * @param inst instance
 * @return 0 if the input data is read successfully, 1 otherwise
 */
int load_instance(instance *inst)
{
	// No need to check anything, we already do it parsing command line arguments
	if (strcmp(inst->input_file, "NULL") == 0)
	{
		generate_random_nodes(inst, inst->nnodes, inst->randomseed);
	}
	else
	{

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

		int active_section = 0; // =0 GENERAL SETTINGS, =1 NODE_COORD_SECTION

		// int do_print = ( VERBOSE >= 1000 );

		while (fgets(line, sizeof(line), fin) != NULL)
		{
			if (VERBOSE >= 2000)
			{
				printf("%s", line);
				fflush(NULL);
			}
			if (strlen(line) <= 1)
				continue; // skip empty lines
			par_name = strtok(line, " :");
			if (VERBOSE >= 3000)
			{
				printf("parameter \"%s\" ", par_name);
				fflush(NULL);
			}

			if (strncmp(par_name, "NAME", 4) == 0)
			{
				active_section = 0;
				continue;
			}

			if (strncmp(par_name, "COMMENT", 7) == 0)
			{
				active_section = 0;
				token1 = strtok(NULL, "");
				// if ( VERBOSE >= 10 ) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
				continue;
			}

			if (strncmp(par_name, "TYPE", 4) == 0)
			{
				token1 = strtok(NULL, " :");
				if (strncmp(token1, "TSP", 3) != 0)
				{
					print_error(" format error:  only TYPE == TSP implemented so far!!!!!!");
					return EXIT_FAILURE;
				}
				active_section = 0;
				continue;
			}

			if (strncmp(par_name, "DIMENSION", 9) == 0)
			{
				if (inst->nnodes >= 0)
				{
					print_error(" repeated DIMENSION section in input file");
					return EXIT_FAILURE;
				}
				token1 = strtok(NULL, " :");
				inst->nnodes = atoi(token1);

				inst->xcoord = (double *)calloc(inst->nnodes, sizeof(double));
				inst->ycoord = (double *)calloc(inst->nnodes, sizeof(double));
				inst->cost_matrix = (double **)calloc(inst->nnodes, sizeof(double));
				for (int i = 0; i < inst->nnodes; i++)
				{
					inst->cost_matrix[i] = (double *)calloc(inst->nnodes, sizeof(double));
				}

				active_section = 0;
				continue;
			}

			if (strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0)
			{
				token1 = strtok(NULL, " :");
				if (strncmp(token1, "ATT", 3) != 0)
				{
					print_error(" format error:  only EDGE_WEIGHT_TYPE == ATT implemented so far!!!!!!");
					return EXIT_FAILURE;
				}
				active_section = 0;
				continue;
			}

			if (strncmp(par_name, "NODE_COORD_SECTION", 18) == 0)
			{
				if (inst->nnodes <= 0)
				{
					print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
					return EXIT_FAILURE;
				}
				active_section = 1;
				continue;
			}

			if (strncmp(par_name, "EOF", 3) == 0)
			{
				active_section = 0;
				break;
			}

			if (active_section == 1) // within NODE_COORD_SECTION
			{
				int i = atoi(par_name) - 1;
				if (i < 0 || i >= inst->nnodes)
				{
					print_error(" ... unknown node in NODE_COORD_SECTION section");
					return EXIT_FAILURE;
				}
				token1 = strtok(NULL, " :,");
				token2 = strtok(NULL, " :,");
				inst->xcoord[i] = atof(token1);
				inst->ycoord[i] = atof(token2);
				// if ( do_print ) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i+1, inst->xcoord[i], inst->ycoord[i]);
				continue;
			}

			return EXIT_FAILURE;
		}

		fclose(fin);
	}

	// Compute cost matrix
	if (tsp_compute_costs(inst) != 0)
	{
		print_error("Error computing cost matrix");
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

/**
 * Print error message
 * @param err_message error message
 * @return void
 */
void print_error(const char *err_message) { fprintf(stderr, "ERROR: %s \n\n", err_message); }

/**
 * Parse command line
 * @param argc number of arguments
 * @param argv arguments
 * @param inst instance
 * @return 0 if the command line is parsed successfully, 1 otherwise
 */
int parse_command_line(int argc, char **argv, instance *inst)
{
	if (VERBOSE >= 100)
		printf("Running %s with %d parameters\n", argv[0], argc - 1);

	int help = 0;
	int source = -1; // 0 for rand, 1 for input file
	int isset_seed = 0;
	int isset_nnodes = 0;

	if (argc < 1)
	{
		help = 1;
	}

	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "--random") || !strcmp(argv[i], "-r"))
		{
			if (source != -1)
			{
				fprintf(stderr, "Error: Cannot mix --random and --file modes\n");
				help = 1;
			}
			source = 0;
		}
		else if (!strcmp(argv[i], "--seed") || !strcmp(argv[i], "-s"))
		{
			if (i + 1 < argc)
			{
				inst->randomseed = atoi(argv[++i]);
				if (inst->randomseed <= 0)
				{
					fprintf(stderr, "Error: Seed must be a positive integer\n");
					help = 1;
				}
				isset_seed = 1;
			}
			else
			{
				fprintf(stderr, "Error: Seed value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--nnodes") || !strcmp(argv[i], "-n"))
		{
			if (i + 1 < argc)
			{
				inst->nnodes = atoi(argv[++i]);
				if (inst->nnodes <= 0)
				{
					fprintf(stderr, "Error: Number of nodes must be a positive integer\n");
					help = 1;
				}
				isset_nnodes = 1;
			}
			else
			{
				fprintf(stderr, "Error: Number of nodes is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--file") || !strcmp(argv[i], "-f"))
		{
			if (source != -1)
			{
				fprintf(stderr, "Error: Cannot mix --random and --file modes\n");
				help = 1;
			}

			if (i + 1 < argc)
			{
				snprintf(inst->input_file, sizeof(inst->input_file), "storage/%s", argv[++i]);
				source = 1;
			}
			else
			{
				fprintf(stderr, "Error: Input file name is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--time_limit") || !strcmp(argv[i], "-tl"))
		{
			if (i + 1 < argc)
			{
				inst->timelimit = atof(argv[++i]);
				if (inst->timelimit <= 0)
				{
					fprintf(stderr, "Error: Time limit must be a positive number\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Time limit value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--method") || !strcmp(argv[i], "-m"))
		{
			if (i + 1 < argc)
			{
				snprintf(inst->method, sizeof(inst->method), "%s", argv[++i]);
				inst->nmethods = 1; // At least one method
				for (char *c = inst->method; *c != '\0'; c++)
				{
					if (*c == '+')
					{
						inst->nmethods++;
					}
				}
			}
			else
			{
				fprintf(stderr, "Error: Method value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--plot"))
		{
			if (i < argc)
			{
				inst->plot = 0;
				i++;
			}
			else
			{
				perror("Plot is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
		{
			help = 1;
			break;
		}
		else if (!strcmp(argv[i], "--vns_kmin") || !strcmp(argv[i], "-vkmin"))
		{
			if (i + 1 < argc)
			{
				inst->vns_kmin = atoi(argv[++i]);
				if (inst->vns_kmin < 1)
				{
					fprintf(stderr, "vns_kmin must be at least 1\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: VNS k min value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--vns_kmax") || !strcmp(argv[i], "-vkmax"))
		{
			if (i + 1 < argc)
			{
				inst->vns_kmax = atoi(argv[++i]);
				if (inst->vns_kmax < 1)
				{
					fprintf(stderr, "vns_kmax must be at least 1\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: VNS k max value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--vns_lr") || !strcmp(argv[i], "-vlr"))
		{
			if (i + 1 < argc)
			{
				inst->vns_learning_rate = atof(argv[++i]);
				if (inst->vns_learning_rate <= 0)
				{
					fprintf(stderr, "Error: VNS learning rate must be positive\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: VNS learning rate value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--vns_jumps") || !strcmp(argv[i], "-vjps"))
		{
			if (i + 1 < argc)
			{
				inst->vns_jumps = atoi(argv[++i]);
				if (inst->vns_jumps < 0)
				{
					fprintf(stderr, "Error: VNS jumps must be non-negative\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: VNS jumps value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--tabu-tenure") || !strcmp(argv[i], "-tten"))
		{
			if (i + 1 < argc)
			{
				inst->tabu_tenure = atoi(argv[++i]);

				if (inst->tabu_tenure < 0)
				{
					fprintf(stderr, "Error: TABU initial tenure size must be non-negative\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Tabu tenure value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--tabu-min") || !strcmp(argv[i], "-ttenmin"))
		{
			if (i + 1 < argc)
			{
				inst->tabu_min = atoi(argv[++i]);

				if (inst->tabu_min < 0)
				{
					fprintf(stderr, "Error: TABU min tenure size must be non-negative\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Tabu min tenure value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--tabu-max") || !strcmp(argv[i], "-ttenmax"))
		{
			if (i + 1 < argc)
			{
				inst->tabu_max = atoi(argv[++i]);

				if (inst->tabu_max < 0)
				{
					fprintf(stderr, "Error: TABU max tenure size must be non-negative\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Tabu max tenure value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--tabu-noimprove") || !strcmp(argv[i], "-tnoimpr"))
		{
			if (i + 1 < argc)
			{
				inst->tabu_noimprove = atoi(argv[++i]);

				if (inst->tabu_noimprove < 0)
				{
					fprintf(stderr, "Error: TABU no improve limit must be non-negative\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Tabu no improve value is missing\n");
				help = 1;
			}
		}
	}

	// Validation
	if (source == 0) // random
	{
		if (!isset_seed || !isset_nnodes)
		{
			fprintf(stderr, "Error: For random generation, both --seed and --nnodes are required\n");
			help = 1;
		}
	}
	else if (source == 1) // file
	{
		if (inst->input_file[0] == '\0')
		{
			fprintf(stderr, "Error: No input file specified\n");
			help = 1;
		}
	}
	else if (source == -1)
	{
		fprintf(stderr, "Error: You must specify either --random or --file mode\n");
		help = 1;
	}

	if (inst->vns_kmin > inst->vns_kmax)
	{
		fprintf(stderr, "Error: vns_kmin must be less than or equal to vns_kmax\n");
		help = 1;
	}

	if (inst->tabu_min > inst->tabu_max || inst->tabu_tenure > inst->tabu_max || inst->tabu_min > inst->tabu_tenure)
	{
		fprintf(stderr, "Error: tabu params are not coherent.\n");
		help = 1;
	}

	if (VERBOSE >= 10)
	{
		printf("\n===== PARAMETERS SET =====\n");
		printf("%-15s : %s\n", "--file", inst->input_file[0] ? inst->input_file : "(none)");
		printf("%-15s : %s\n", "--method", inst->method);
		if (inst->nnodes > 0)
			printf("%-15s : %d\n", "--nnodes", inst->nnodes);
		printf("%-15s : %d\n", "--seed", inst->randomseed);
		printf("%-15s : %.2f\n", "--time_limit", inst->timelimit);
		printf("===========================\n\n");
	}

	if (help)
	{
		printf("Options:\n");
		printf("  -r, --random              Generate random nodes\n");
		printf("  -n, --nnodes <int>        Number of nodes (required with -r)\n");
		printf("  -s, --seed <int>          Random seed (required with -r)\n");
		printf("  -f, --file <filename>     Input file to load instance\n");
		printf("  -m, --method <string>     Method to solve TSP (e.g., n_n, n_n+apply_two_opt)\n");
		printf("  -p, --plot                Paramater to see and save the solution plot\n");

		printf("  -tl,--time_limit <float>  Time limit in seconds\n");
		printf("  -h, --help                Show this help message\n");

		printf("\nExamples:\n");
		printf("  ./main -r -n 50 -s 123 -m n_n -tl 10\n");
		printf("  ./main -f instance.txt -m random -tl 20\n\n");

		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

/**
 * Generate random nodes
 * @param inst instance
 * @param nnodes number of nodes
 * @param seed seed
 * @return void
 */
void generate_random_nodes(instance *inst, int nnodes, int seed)
{
	if (VERBOSE >= 50)
		printf("Generating random nodes...\n");

	for (int i = 0; i < nnodes; i++)
	{
		inst->xcoord[i] = rand01() * MAX_X;
		inst->ycoord[i] = rand01() * MAX_Y;
	}

	if (VERBOSE >= 50)
		printf("Random nodes generated\n");
}

/**
 * Free instance
 * @param inst instance
 * @return void
 */
void free_instance(instance *inst)
{
	free(inst->xcoord);
	inst->xcoord = NULL;

	free(inst->ycoord);
	inst->ycoord = NULL;

	for (int i = 0; i < inst->nnodes; i++)
	{
		free(inst->cost_matrix[i]);
	}
	free(inst->cost_matrix);
	inst->cost_matrix = NULL;
}

/**
 * Free solution
 * @param sol solution
 * @return void
 */
void free_sol(solution *sol)
{
	if (sol && sol->tour)
	{
		free(sol->tour);
		sol->tour = NULL;
	}
}

/**
 * Generate a random path using the Fisher-Yates shuffle algorithm
 * @param inst the instance
 * @param sol path of the solution
 * @return 0 if the random path is generated successfully, 1 otherwise
 */
int generate_random_path(const instance *inst, solution *sol)
{
	int nnodes = inst->nnodes;
	int nodes[nnodes];
	for (int i = 0; i < nnodes; i++)
	{
		nodes[i] = i;
	}

	for (int i = nnodes - 1; i > 0; i--)
	{
		int j = rand() % (i + 1);
		// Swap nodes[i] with nodes[j]
		int temp = nodes[i];
		nodes[i] = nodes[j];
		nodes[j] = temp;
	}

	for (int i = 0; i < nnodes; i++)
	{
		sol->tour[i] = nodes[i];
	}
	sol->tour[nnodes] = sol->tour[0]; // Comes back to the initial node

	evaluate_path_cost(inst, sol);

	return EXIT_SUCCESS;
}

/**
 * Print path
 * @param inst instance
 * @param sol solution path
 * @return void
 */
void print_solution_path(const instance *inst, const solution *sol)
{
	for (int i = 0; i <= inst->nnodes; i++)
	{
		int node = sol->tour[i];
		printf("(%d)", node);
	}
	printf("\n");
}

/**
 * Print nodes
 * @param inst instance
 * @return void
 */
void print_node_coordinates(instance *inst)
{
	for (int i = 0; i < inst->nnodes; i++)
	{
		printf("Node %d: (%d, %d)\n", i, (int)inst->xcoord[i], (int)inst->ycoord[i]);
	}
}

/**
 * Check time
 * @param inst instance
 * @param time starting_time
 * @return 0 if the time limit is not reached, 1 otherwise
 */
int check_time(const instance *inst, double starting_time)
{
	double current_time = second();
	if (current_time - starting_time > inst->timelimit)
	{
		// print_error("Time limit reached\n");
		return EXIT_FAILURE;
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
int write_path_to_file(const instance *inst, const solution *sol, const char *filename)
{
	FILE *fout = fopen(filename, "w");
	if (fout == NULL)
	{
		print_error("File not found");
		return EXIT_FAILURE;
	}

	for (int i = 0; i <= inst->nnodes; i++)
	{
		int node = sol->tour[i];
		fprintf(fout, "%d %d\n", (int)inst->xcoord[node], (int)inst->ycoord[node]);
	}

	fclose(fout);
	return EXIT_SUCCESS;
}

/**
 * cost of the path
 * @param inst instance
 * @param sol solution path
 * @return 0 if the cost of the path is calculated successfully, 1 otherwise
 */
int evaluate_path_cost(const instance *inst, solution *sol)
{
	double cost_p = 0;
	for (int i = 0; i < inst->nnodes - 1; i++)
	{
		double edge_cost = inst->cost_matrix[sol->tour[i]][sol->tour[i + 1]];

		if (VERBOSE >= 70)
			printf("Cost from node %d to node %d: %lf\n", sol->tour[i], sol->tour[i + 1], edge_cost);
		cost_p += edge_cost;
	}

	// Add the cost to return to the starting node
	double return_cost = inst->cost_matrix[sol->tour[inst->nnodes - 1]][sol->tour[0]];

	if (VERBOSE >= 70)
		printf("Cost from node %d to node %d: %lf\n", sol->tour[inst->nnodes - 1], sol->tour[0], return_cost);
	cost_p += return_cost;
	// printf("Total cost: %lf\n", cost_p);
	sol->cost = cost_p; // update the cost of the solution

	return EXIT_SUCCESS;
}

/**
 * Two-opt heuristic
 * @param inst instance
 * @param sol solution path
 * @return 0 if the two-opt heuristic is applied successfully, 1 otherwise
 */
int apply_two_opt(const instance *inst, solution *sol)
{
	if (sol == NULL) // Check if sol is NULL
	{
		print_error("sol is NULL in two_opt");
		return EXIT_FAILURE;
	}

	double starting_time = second();

	if (!sol->initialized) // Check if solution is not initialized
	{
		// Generate a random path if the solution is not initialized
		if (generate_random_path(inst, sol) != EXIT_SUCCESS)
		{
			print_error("Failed to generate random path in two_opt");
			return EXIT_FAILURE;
		}
	}

	double time = second();
	if (VERBOSE >= 70)
	{
		printf("Applying two-opt heuristic...\n");
	}
	int nnodes = inst->nnodes;
	int improved = 1;			  // flag to indicate if the path has been improved
	double best_cost = sol->cost; // Initialize with the current cost
	while (improved)
	{
		if (check_time(inst, starting_time))
			return EXIT_SUCCESS;

		improved = 0; // reset the flag
		for (int i = 0; i < nnodes - 1; i++)
		{
			for (int j = i + 1; j < nnodes; j++)
			{
				double delta_cost = path_cost_delta(i, j, sol, inst);
				if (delta_cost < 0)
				{
					reverse_path_segment(i, j, sol); // swap nodes i and j
					improved = 1;					 // set the flag
					sol->cost += delta_cost;		 // Update the cost
					if (sol->cost < best_cost)
					{
						best_cost = sol->cost; // Update the best cost
					}
				}
			}
		}
	}

	if (VERBOSE >= 70)
	{
		printf("Two-opt heuristic applied\n");
		printf("Elapsed time: %lf\n", second() - time);
		printf("--------------------------------------------\n");
	}

	if (!sol->initialized) // Mark as initialized only if not already initialized
	{
		sol->initialized = 1;
	}

	evaluate_path_cost(inst, sol);

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
double path_cost_delta(int i, int j, const solution *sol, const instance *inst)
{
	return (inst->cost_matrix[sol->tour[i + 1]][sol->tour[j + 1]] + inst->cost_matrix[sol->tour[i]][sol->tour[j]] - (inst->cost_matrix[sol->tour[i]][sol->tour[i + 1]] + inst->cost_matrix[sol->tour[j]][sol->tour[j + 1]]));
}

/**
 * Swap nodes i and j in the path
 * @param i node i
 * @param j node j
 * @param sol solution path
 * @return void
 */
void reverse_path_segment(int i, int j, solution *sol)
{
	while (++i < j)
	{
		int temp = sol->tour[i];
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
int execute_selected_method(instance *inst, solution *sol)
{
	// Initialized the sol
	sol->tour = (int *)malloc((inst->nnodes + 1) * sizeof(int));
	if (sol->tour == NULL)
	{
		print_error("Memory allocation failed for sol->tour");
		return EXIT_FAILURE;
	}
	sol->initialized = 0; // Initialized to 0 (false)

	char *method_str = strdup(inst->method); // Double the string to not modify it directly
	if (method_str == NULL)
	{
		print_error("Memory allocation failed");
		free(sol->tour);
		return EXIT_FAILURE;
	}

	char *method = strtok(method_str, "+"); // Divide the string in separate methods by the '+'

	if (method != NULL)
		inst->timelimit = inst->timelimit / inst->nmethods; // Time limit for each method

	while (method != NULL)
	{
		double starting_time_method = second();
		if (VERBOSE >= 50)
		{
			printf("\n----------------------------------\n");
			printf("METHOD %s \n\n", method);
		}

		if (strcmp(method, "n_n") == 0)
		{
			if (apply_greedy_search(inst, sol)) // Nearest neighbor heuristic
			{
				print_error("Error applying nearest neighbor heuristic");
				free(method_str);
				return EXIT_FAILURE;
			}
			if (VERBOSE >= 50)
				printf("Nearest Neighbor done in %lf seconds\n\n", second() - starting_time_method);
		}
		else if (strcmp(method, "two_opt") == 0)
		{
			if (apply_two_opt(inst, sol))
			{
				print_error("Error applying 2-opt");
				free(method_str);
				return EXIT_FAILURE;
			}

			if (VERBOSE >= 50)
				printf("Two Opt done in %lf seconds\n\n", second() - starting_time_method);
		}
		else if (strcmp(method, "random") == 0)
		{
			if (generate_random_path(inst, sol))
			{
				print_error("Error generating random path");
				free(method_str);
				return EXIT_FAILURE;
			}

			if (VERBOSE >= 50)
				printf("Random path done in %lf seconds\n\n", second() - starting_time_method);
		}
		else if (strcmp(method, "vns") == 0)
		{
			if (apply_greedy_search(inst, sol)) // Nearest neighbor heuristic
			{
				print_error("Error applying nearest neighbor heuristic");
				free(method_str);
				return EXIT_FAILURE;
			}
			if (apply_heuristic_vns(inst, sol))
			{
				print_error("Error applying VNS");
				free(method_str);
				return EXIT_FAILURE;
			}

			if (VERBOSE >= 50)
				printf("VNS done in %lf seconds\n\n", second() - starting_time_method);
		}
		else if (strcmp(method, "tabu") == 0)
		{
			if (apply_greedy_search(inst, sol)) // Nearest neighbor heuristic
			{
				print_error("Error applying nearest neighbor heuristic");
				free(method_str);
				return EXIT_FAILURE;
			}
			if (apply_heuristic_tabu(inst, sol))
			{
				print_error("Error applying Tabu Search");
				free(method_str);
				return EXIT_FAILURE;
			}

			if (VERBOSE >= 50)
				printf("Tabu search done in %lf seconds\n\n", second() - starting_time_method);
		}
		else if (strcmp(method, "benders") == 0)
		{
			if (apply_cplex_beneders(inst, sol))
			{
				print_error("Error Solving with CPLEX");
				free(method_str);
				return EXIT_FAILURE;
			}

			if (VERBOSE >= 50)
				printf("CPLEX done in %lf seconds\n\n", second() - starting_time_method);
		}
		else
		{
			print_error("Invalid method");
			print_error("USAGE: -method [n_n|two_opt|random|vns|tabu] (separated by '+')");
			free(method_str);
			return EXIT_FAILURE;
		}

		/*
		if (check_time(inst, starting_time))
		{
			free(method_str);
			evaluate_path_cost(inst, sol);
			return EXIT_SUCCESS;
		}
			*/

		evaluate_path_cost(inst, sol);
		method = strtok(NULL, "+"); // Get the next method
	}

	free(method_str);

	// Save the cost of the best solution
	evaluate_path_cost(inst, sol);
	return EXIT_SUCCESS;
}

void allocate_buffers(instance *tsp)
{
	// Libera se già allocato (difensivo, opzionale se sicuro altrove)
	if (tsp->xcoord)
	{
		free(tsp->xcoord);
	}
	tsp->xcoord = (double *)calloc(tsp->nnodes, sizeof(double));
	if (!tsp->xcoord)
		print_error("Failed to allocate xcoord");

	if (tsp->ycoord)
	{
		free(tsp->ycoord);
	}
	tsp->ycoord = (double *)calloc(tsp->nnodes, sizeof(double));
	if (!tsp->ycoord)
		print_error("Failed to allocate ycoord");

	if (tsp->cost_matrix)
	{
		for (int i = 0; i < tsp->nnodes; i++)
		{
			if (tsp->cost_matrix[i])
			{
				free(tsp->cost_matrix[i]);
			}
		}
		free(tsp->cost_matrix);
	}
	tsp->cost_matrix = (double **)calloc(tsp->nnodes, sizeof(double *));
	if (!tsp->cost_matrix)
		print_error("Failed to allocate cost_matrix");

	for (int i = 0; i < tsp->nnodes; i++)
	{
		tsp->cost_matrix[i] = (double *)calloc(tsp->nnodes, sizeof(double));
		if (!tsp->cost_matrix[i])
			print_error("Failed to allocate cost_matrix row");
	}
}

void init(instance *inst)
{
	memset(inst, 0, sizeof(instance));

	inst->nnodes = 0;
	inst->xcoord = NULL;
	inst->ycoord = NULL;

	inst->cost_matrix = NULL;
	// inst->best_sol.tour = NULL;
	// inst->best_sol.cost = CPX_INFBOUND;

	strcpy(inst->input_file, "NULL");
	strcpy(inst->method, "NULL");
	inst->timelimit = CPX_INFBOUND;
	inst->randomseed = rand();
	inst->plot = 1;

	inst->starting_time = second();

	inst->vns_kmin = 1;
	inst->vns_kmax = 5;
	inst->vns_learning_rate = 1.0;
	inst->vns_jumps = 0;

	inst->tabu_min = 0;
	inst->tabu_max = 0;
	inst->tabu_tenure = 0;
	inst->tabu_noimprove = 0;

	strcpy(inst->csv_filename, "");
}

/**
 * Generate 3-opt positions: fills positions[0..2] with valid indices for a 3-opt move.
 * @param tsp TSP instance
 * @param positions array of 3 integers
 * @return void
 */
static void generate_three_opt_positions(instance *tsp, int *positions)
{
	while (1)
	{
		positions[0] = rand() % tsp->nnodes;
		positions[1] = rand() % tsp->nnodes;
		positions[2] = rand() % tsp->nnodes;
		qsort(positions, 3, sizeof(int), compar);

		// Ensure the positions are sufficiently apart:
		if (positions[1] > positions[0] + 1 &&
			positions[2] > positions[1] + 1 &&
			positions[2] + 1 < tsp->nnodes)
		{
			break;
		}
	}
}

/**
 * Perform a 3-opt swap on the current solution, storing the result in the new solution.
 * @param tsp TSP instance
 * @param current_sol current solution
 * @param new_tour new solution path
 * @param positions array of 3 integers
 * @return void
 */
static void perform_three_opt_swap(instance *tsp, solution *current_sol, int *new_tour, int *positions)
{
	// Copy segment from index 0 to positions[0] unchanged.
	memcpy(new_tour, current_sol->tour, sizeof(int) * (positions[0] + 1));
	int pos = positions[0] + 1;

	// Reverse the segment between positions[0]+1 and positions[1].
	for (int i = positions[1]; i >= positions[0] + 1; i--)
	{
		new_tour[pos++] = current_sol->tour[i];
	}

	// Reverse the segment between positions[1]+1 and positions[2].
	for (int i = positions[2]; i >= positions[1] + 1; i--)
	{
		new_tour[pos++] = current_sol->tour[i];
	}

	// Copy the remainder of the tour from positions[2]+1 to the end.
	memcpy(new_tour + pos, current_sol->tour + positions[2] + 1, sizeof(int) * (tsp->nnodes - positions[2] - 1));
}

/**
 * Recompute the solution cost after applying a 3-opt swap.
 * @param tsp TSP instance
 * @param solution solution path
 * @return 0 if the solution cost is recomputed successfully, -1 otherwise
 */
int apply_three_opt(instance *tsp, solution *sol)
{
	if (!tsp->cost_matrix || tsp->nnodes <= 0)
	{
		return -1;
	}

	int positions[3];
	generate_three_opt_positions(tsp, positions);

	int *temp_tour = (int *)malloc(sizeof(int) * (tsp->nnodes + 1));
	if (temp_tour == NULL)
	{
		return -1;
	}

	// Apply the 3-opt swap using the generated positions.
	perform_three_opt_swap(tsp, sol, temp_tour, positions);
	memcpy(sol->tour, temp_tour, sizeof(int) * (tsp->nnodes + 1));

	free(temp_tour);

	// Recompute and update the solution cost.
	evaluate_path_cost(tsp, sol);

	return 0;
}
