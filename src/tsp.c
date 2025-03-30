#include "tsp.h"
#include "utils.h"
#include "tsp_heuristics.h"
#include "tsp_greedy.h"

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
			}
			else
			{
				fprintf(stderr, "Error: Method value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
		{
			help = 1;
			break;
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

	if (VERBOSE >= 10)
	{
		printf("\n===== PARAMETERS SET =====\n");
		printf("%-15s : %s\n", "--file", inst->input_file[0] ? inst->input_file : "(none)");
		printf("%-15s : %s\n", "--method", inst->method);
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

	set_seed(seed);
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
void free_instance(instance *inst, solution *sol)
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

	if (sol && sol->tour)
	{
		free(sol->tour);
		sol->tour = NULL;
	}
}

/**
 * Generate a random path using the Fisher-Yates shuffle algorithm
 * @param sol path of the solution
 * @param nnodes number of nodes
 * @param seed seed
 * @return 0 if the random path is generated successfully, 1 otherwise
 */
int generate_random_path(solution *sol, int nnodes, int seed)
{
	int nodes[nnodes];
	for (int i = 0; i < nnodes; i++)
	{
		nodes[i] = i;
	}

	set_seed(seed);
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
 * @param start_time start time
 * @return 0 if the time limit is not reached, 1 otherwise
 */
int check_time(const instance *inst)
{
	double current_time = second();
	if (current_time - inst->starting_time > inst->timelimit)
	{
		print_error("Time limit reached");
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
	double time = second();
	if (VERBOSE >= 50)
	{
		printf("Applying two-opt heuristic...\n");
	}
	int nnodes = inst->nnodes;
	int improved = 1; // flag to indicate if the path has been improved
	while (improved)
	{
		improved = 0; // reset the flag
		for (int i = 0; i < nnodes - 1; i++)
		{
			for (int j = i + 1; j < nnodes; j++)
			{
				if (path_cost_delta(i, j, sol, inst) < 0)
				{									 // if the path is improved
					reverse_path_segment(i, j, sol); // swap nodes i and j
					improved = 1;					 // set the flag
				}
			}
		}
	}

	if (VERBOSE >= 50)
	{
		printf("Two-opt heuristic applied\n");
		printf("Elapsed time: %lf\n", second() - time);
		printf("--------------------------------------------\n");
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
	if (strcmp(inst->method, "n_n") == 0)
	{
		if (apply_greedy_search(inst, sol)) // Nearest neighbor heuristic
		{
			print_error("Error applying nearest neighbor heuristic");
			return EXIT_FAILURE;
		}
		if (check_time(inst))
		{
			return EXIT_FAILURE;
		}
	}
	else if (strcmp(inst->method, "n_n+two_opt") == 0)
	{
		if (solve_greedy(inst, sol)) // Nearest neighbor heuristic with two_opt
		{
			print_error("Nearest neighbor heuristic failed");
			return EXIT_FAILURE;
		}
		if (apply_two_opt(inst, sol))
		{
			print_error("Error applying 2-opt");
			return EXIT_FAILURE;
		}
		if (check_time(inst))
		{
			return EXIT_FAILURE;
		}
	}
	else if (strcmp(inst->method, "random+two_opt") == 0)
	{ // Random path with two_opt
		if (random_path(sol, inst->nnodes, inst->randomseed))
		{
			print_error("Error generating random path");
			return EXIT_FAILURE;
		}
		if (apply_two_opt(inst, sol))
		{
			print_error("Error applying 2-opt");
			return EXIT_FAILURE;
		}
		if (check_time(inst))
		{
			return EXIT_FAILURE;
		}
	}
	else if (strcmp(inst->method, "random") == 0)
	{
		generate_random_path(sol, inst->nnodes, inst->randomseed); // Random path
		if (check_time(inst))
		{
			print_error("VNS failed");
			return EXIT_FAILURE;
		}
	}
	else if (strcmp(inst->method, "vns") == 0)
	{ // VNS
		if (apply_heuristic_vns(inst, sol))
		{
			print_error("Error applying VNS");
			return EXIT_FAILURE;
		}
	}
	else if (strcmp(inst->method, "tabu") == 0) // TABU SEARCH
	{
		if (apply_heuristic_tabu(inst, sol))
		{
			print_error("Error applying Tabu Search");
			return EXIT_FAILURE;
		}
	}
	else
	{
		print_error("Invalid method\n");
		print_error("USAGE: -method [n_n|n_n+two_opt|random+two_opt|random]");
		return EXIT_FAILURE;
	}

	// Print the cost of the best solution
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
	inst->randomseed = 0;

	inst->starting_time = second();
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
