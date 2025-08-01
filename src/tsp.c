#include "tsp.h"
#include "utils.h"
#include "tsp_heuristics.h"
#include "tsp_greedy.h"
#include "tsp_cplex.h"
#include <string.h>

/**
 * Load TSP instance data, either from a TSPLIB-style file or by generating
 * random coordinates.
 *
 * @param inst instance
 * @return 0 if the input data is read successfully, 1 otherwise
 */
int load_instance(instance *inst)
{
	if (strcmp(inst->input_file, "NULL") == 0)
	{
		generate_random_nodes(inst, inst->nnodes);
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

		int active_section = 0; // = 0 GENERAL SETTINGS, = 1 NODE_COORD_SECTION

		while (fgets(line, sizeof(line), fin) != NULL)
		{
			if (VERBOSE >= 2000)
			{
				printf("%s", line);
				fflush(NULL);
			}

			if (strlen(line) <= 1)
			{
				continue;
			}

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
				if (strncmp(token1, "ATT", 3) == 0)
					inst->edge_weight = 0;
				else if (strncmp(token1, "EUC_2D", 6) == 0)
					inst->edge_weight = 1;
				else
					inst->edge_weight = -1;

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

			if (active_section == 1)
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

				continue;
			}

			return EXIT_FAILURE;
		}

		fclose(fin);
	}

	if (tsp_compute_costs(inst) != 0)
	{
		print_error("Error computing cost matrix");
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

/**
 * Print an error message to stderr, prefixed with "ERROR:".
 *
 * @param err_message error message
 * @return void
 */
void print_error(const char *err_message) { fprintf(stderr, "ERROR: %s \n\n", err_message); }

/**
 * Parse command-line flags and populate the instance with the data provided.
 *
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
	int source = -1; // = 0 rand, = 1 input file
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
		else if (!strcmp(argv[i], "--greedy-starts") || !strcmp(argv[i], "-gs"))
		{
			if (i + 1 < argc)
			{
				inst->greedy_starts = atoi(argv[++i]);

				if (inst->greedy_starts < 0)
				{
					fprintf(stderr, "Error: Greedy Starts value must be non-negative\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Greedy Starts value is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--hard_fixing_percentage") || !strcmp(argv[i], "-hf_p"))
		{
			if (i + 1 < argc)
			{
				inst->hard_fixing_percentage = atof(argv[++i]);

				if (inst->hard_fixing_percentage < 0 || inst->hard_fixing_percentage > 1)
				{
					fprintf(stderr, "Error: Hard fixing percentage must be between 0 and 1\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Hard fixing percentage is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--hard_fixing_local_time") || !strcmp(argv[i], "-hf_lt"))
		{
			if (i + 1 < argc)
			{
				inst->hard_fixing_local_time = atof(argv[++i]);

				if (inst->hard_fixing_percentage < 0)
				{
					fprintf(stderr, "Error: Hard fixing local time limit must be positive\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Hard fixing local time limit is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--local_branch_k") || !strcmp(argv[i], "-lb_k"))
		{
			if (i + 1 < argc)
			{
				inst->local_branch_k = atof(argv[++i]);

				if (inst->local_branch_k < 0)
				{
					fprintf(stderr, "Error: Local Branching initial K must be positive\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Local Branching initial K is missing\n");
				help = 1;
			}
		}
		else if (!strcmp(argv[i], "--local_branch_step") || !strcmp(argv[i], "-lb_s"))
		{
			if (i + 1 < argc)
			{
				inst->local_branch_step = atof(argv[++i]);

				if (inst->local_branch_step < 0)
				{
					fprintf(stderr, "Error: Local Branching K step must be positive\n");
					help = 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: Local Branching K step is missing\n");
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
		if (strstr(inst->method, "n_n"))
			printf("%-15s : %d\n", "--greedy_starts", inst->greedy_starts);
		if (strstr(inst->method, "vns"))
		{
			printf("%-15s : %d\n", "--vns_kmin", inst->vns_kmin);
			printf("%-15s : %d\n", "--vns_kmax", inst->vns_kmax);
			printf("%-15s : %.2f\n", "--vns_lr", inst->vns_learning_rate);
			printf("%-15s : %d\n", "--vns_jumps", inst->vns_jumps);
		}
		if (strstr(inst->method, "tabu"))
		{
			printf("%-15s : %d\n", "--tabu-tenure", inst->tabu_tenure);
			printf("%-15s : %d\n", "--tabu-min", inst->tabu_min);
			printf("%-15s : %d\n", "--tabu-max", inst->tabu_max);
			printf("%-15s: %d\n", "--tabu-noimprove", inst->tabu_noimprove);
		}
		if (strstr(inst->method, "hard_fix"))
		{
			printf("%-15s : %.2f\n", "--hard_fixing_percentage", inst->hard_fixing_percentage);
			printf("%-15s : %.2f\n", "--hard_fixing_local_time", inst->hard_fixing_local_time);
		}
		if (strstr(inst->method, "local_branch"))
		{
			printf("%-15s : %d\n", "--local_branch_k", inst->local_branch_k);
			printf("%-15s : %d\n", "--local_branch_step", inst->local_branch_step);
		}
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
		printf("  -rs, --runs <int> 	    Number of runs (default: 1)\n");
		printf("  -h, --help                Show this help message\n");

		printf("\nMethods parameters:\n");
		printf("Greedy:\n");
		printf("  -gs, 	  --greedy_starts <int> Number of greedy starts (default: 10)\n");
		printf("VNS:\n");
		printf("  -vkmin, --vns_kmin <int>   Minimum k for VNS (default: 1)\n");
		printf("  -vkmax, --vns_kmax <int>   Maximum k (default: 5)\n");
		printf("  -vlr,   --vns_lr <float>   Learning rate for VNS (default: 1.0)\n");
		printf("  -vjps,  --vns_jumps <int>  Number of jumps for VNS (default: 0)\n");
		printf("Tabu Search:\n");
		printf("  -tten,    --tabu-tenure <int> 	Initial tabu tenure size (default: (min_tenure + max_tenure) / 2)\n");
		printf("  -ttenmin, --tabu-min <int> 		Minimum tabu tenure size (default: nnodes / 10)\n");
		printf("  -ttenmax, --tabu-max <int> 		Maximum tabu tenure size (default: nnodes + nnodes / 10)\n");
		printf("  -tnoimpr, --tabu-noimprove <int> 	No improve limit for Tabu (default: 50)\n");
		printf("Hard Fixing CPLEX:\n");
		printf("  -hf_p, --hard_fixing_percentage <float> Hard fixing percentage (default: 0.3)\n");
		printf("  -hf_lt, --hard_fixing_local_time <float> Hard fixing local time limit at each iteration (default: 15sec)\n");
		printf("Local Branching CPLEX:\n");
		printf("  -lb_k, --local_branch_k    <int> Local Branching initial K (default: 5)\n");
		printf("  -lb_s, --local_branch_step <int> Local Branching K step (default: 5)\n");

		printf("\nExamples:\n");
		printf("  ./tsp_solver -r -n 50 -s 123 -m n_n -tl 10\n");
		printf("  ./tsp_solver -f instance.txt -m random -tl 20\n\n");

		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

/**
 * Generate a set of random nodes in a grid MAX_X * MAX_Y
 *
 * @param inst instance
 * @param nnodes number of nodes
 * @return void
 */
void generate_random_nodes(instance *inst, int nnodes)
{
	if (VERBOSE >= 50)
		printf("Generating random nodes...\n");

	inst->xcoord = (double *)malloc(nnodes * sizeof(double));
	inst->ycoord = (double *)malloc(nnodes * sizeof(double));
	inst->cost_matrix = (double **)calloc(inst->nnodes, sizeof(double));
	for (int i = 0; i < inst->nnodes; i++)
	{
		inst->cost_matrix[i] = (double *)calloc(inst->nnodes, sizeof(double));
	}
	if (inst->xcoord == NULL || inst->ycoord == NULL || inst->cost_matrix == NULL)
	{
		print_error("Memory allocation failed");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < nnodes; i++)
	{
		inst->xcoord[i] = rand01() * MAX_X;
		inst->ycoord[i] = rand01() * MAX_Y;
	}

	if (VERBOSE >= 50)
		printf("Random nodes generated\n");
}

/**
 * Free instance data structures
 *
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
 * Free solution data structures
 *
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
 *
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
 * Print the solution path as a sequence of nodes
 *
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
 * Print all nodes with their coordinates
 *
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
 * Check if the execution time limit has been reached
 *
 * @param inst instance
 * @param time starting_time
 * @return 0 if the time limit is not reached, 1 otherwise
 */
int check_time(const instance *inst, double starting_time)
{
	double current_time = second();
	if (current_time - starting_time > inst->timelimit)
	{
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/**
 * Write solution path to file
 *
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
 * Compute the cost of the provided path
 *
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

		if (VERBOSE >= 80)
			printf("Cost from node %d to node %d: %lf\n", sol->tour[i], sol->tour[i + 1], edge_cost);
		cost_p += edge_cost;
	}

	// Add the cost to return to the starting node
	double return_cost = inst->cost_matrix[sol->tour[inst->nnodes - 1]][sol->tour[0]];

	if (VERBOSE >= 80)
		printf("Cost from node %d to node %d: %lf\n", sol->tour[inst->nnodes - 1], sol->tour[0], return_cost);
	cost_p += return_cost;

	sol->cost = cost_p; // update the cost of the solution

	return EXIT_SUCCESS;
}

/**
 * Apply the 2-opt local optimization heuristic to the current tour.
 *
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

	if (!sol->initialized)
	{
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
	int improved = 1;
	double best_cost = sol->cost;
	while (improved)
	{
		if (check_time(inst, starting_time))
			return EXIT_SUCCESS;

		improved = 0;
		for (int i = 0; i < nnodes - 1; i++)
		{
			for (int j = i + 1; j < nnodes; j++)
			{
				double delta_cost = path_cost_delta(i, j, sol, inst);
				if (delta_cost < 0)
				{
					reverse_path_segment(i, j, sol);
					improved = 1;
					sol->cost += delta_cost;
					if (sol->cost < best_cost)
					{
						best_cost = sol->cost;
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

	if (!sol->initialized)
	{
		sol->initialized = 1;
	}

	evaluate_path_cost(inst, sol);

	return EXIT_SUCCESS;
}

/**
 * Calculate the delta in the cost of the path when swapping nodes i and j in the tour
 *
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
 *
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
 *
 * @param tsp instance
 * @return 0 if the cost matrix is computed successfully, 1 otherwise
 */
int tsp_compute_costs(instance *tsp)
{
	if (tsp->cost_matrix == NULL || tsp->xcoord == NULL || tsp->ycoord == NULL)
	{
		print_error("Cost matrix, xcoord or ycoord is NULL");
		return EXIT_FAILURE;
	}

	for (int i = 0; i < tsp->nnodes; i++)
	{
		for (int j = 0; j < tsp->nnodes; j++)
		{
			if (i == j)
			{
				tsp->cost_matrix[i][j] = 0; // Distance to itsself is 0
				continue;
			}

			double deltax = tsp->xcoord[i] - tsp->xcoord[j];
			double deltay = tsp->ycoord[i] - tsp->ycoord[j];
			double dist_squared = deltax * deltax + deltay * deltay;
			double dist;

			if (tsp->edge_weight == 0) // ATT
			{
				dist = sqrt(dist_squared / 10.0);
				int round_d = (int)(dist + 0.5);
				if ((double)round_d < dist)
					round_d++;

				tsp->cost_matrix[i][j] = (int)round_d;
			}
			else if (tsp->edge_weight == 1) // EUC_2D
			{
				dist = sqrt(dist_squared);
				tsp->cost_matrix[i][j] = (int)(dist + 0.5);
			}

			else
			{
				print_error("Unknown edge weight type");
				return EXIT_FAILURE;
			}
		}
	}

	return EXIT_SUCCESS;
}

/**
 * Run the selected set of methods
 *
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
	sol->initialized = 0;

	char *method_str = strdup(inst->method);
	if (method_str == NULL)
	{
		print_error("Memory allocation failed");
		free(sol->tour);
		return EXIT_FAILURE;
	}

	char *method = strtok(method_str, "+"); // Divide the string in separate methods by the '+'

	double total_weight = 0.0;
	double old_timelimit = inst->timelimit;
	double remaining_time = old_timelimit;
	double remaining_weight = 0.0;

	while (method != NULL)
	{
		total_weight += get_method_weight(method);
		method = strtok(NULL, "+");
	}
	free(method_str);

	remaining_weight = total_weight;

	method_str = strdup(inst->method);
	if (method_str == NULL)
	{
		print_error("Memory allocation failed");
		free(sol->tour);
		return EXIT_FAILURE;
	}
	method = strtok(method_str, "+");

	while (method != NULL)
	{
		double weight = get_method_weight(method);
		double slice = remaining_time * (weight / remaining_weight);
		inst->timelimit = slice;

		if (VERBOSE >= 50)
		{
			printf("\n----------------------------------\n");
			printf("METHOD %s \n\n", method);
		}

		printf("Time limit for this method: %.2f s\n", inst->timelimit);

		double starting_time_method = second();

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
		else if (strcmp(method, "two_opt") == 0) // Two-opt heuristic
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
		else if (strcmp(method, "random") == 0) // Random path generation
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
		else if (strcmp(method, "vns") == 0) // Variable Neighbourhood Search (VNS)
		{
			if (apply_greedy_search(inst, sol))
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
		else if (strcmp(method, "tabu") == 0) // Tabu Search
		{
			if (apply_greedy_search(inst, sol))
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
		else if (strcmp(method, "benders") == 0) // Benders Decomposition
		{
			if (apply_cplex_benders(inst, sol))
			{
				print_error("Error Solving with CPLEX");
				free(method_str);
				return EXIT_FAILURE;
			}

			if (VERBOSE >= 50)
				printf("CPLEX done in %lf seconds\n\n", second() - starting_time_method);
		}
		else if (strcmp(method, "branch_and_cut") == 0) // Branch and Cut
		{
			if (apply_cplex_branchcut(inst, sol))
			{
				print_error("Error Solving with CPLEX");
				free(method_str);
				return EXIT_FAILURE;
			}

			if (VERBOSE >= 50)
				printf("CPLEX done in %lf seconds\n\n", second() - starting_time_method);
		}
		else if (strcmp(method, "hard_fix") == 0) // Hard Fixing
		{
			if (apply_cplex_hardfix(inst, sol))
			{
				print_error("Error Solving with CPLEX");
				free(method_str);
				return EXIT_FAILURE;
			}

			if (VERBOSE >= 50)
				printf("CPLEX done in %lf seconds\n\n", second() - starting_time_method);
		}
		else if (strcmp(method, "local_branch") == 0) // Local Branching
		{
			if (apply_cplex_localbranch(inst, sol))
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
			print_error("USAGE: -method [n_n|two_opt|random|vns|tabu|benders|branch_and_cut|hard_fix|local_branch] (separated by '+')");
			free(method_str);
			return EXIT_FAILURE;
		}

		double used = second() - starting_time_method;
		remaining_time -= used;
		remaining_weight -= weight;
		if (remaining_time <= 0.0)
		{
			if (VERBOSE >= 50)
				printf("Time pool exhausted—halting remaining methods\n");
			break;
		}

		evaluate_path_cost(inst, sol);
		method = strtok(NULL, "+");
	}

	free(method_str);

	inst->timelimit = old_timelimit;

	evaluate_path_cost(inst, sol);

	return EXIT_SUCCESS;
}

/**
 * Allocate or reallocate the coordinate and cost‐matrix buffers
 * for an instance, freeing any previous allocations.
 *
 * @param tsp  Pointer to the TSP instance whose buffers should be (re)allocated.
 */
void allocate_buffers(instance *tsp)
{
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

/**
 * Initialize all fields of a TSP instance to default values.
 *
 * @param inst  Pointer to the instance to initialize.
 */
void init(instance *inst)
{
	memset(inst, 0, sizeof(instance));

	inst->nnodes = 0;
	inst->xcoord = NULL;
	inst->ycoord = NULL;

	inst->cost_matrix = NULL;

	strcpy(inst->input_file, "NULL");
	strcpy(inst->method, "NULL");
	inst->timelimit = CPX_INFBOUND;
	inst->randomseed = rand();
	inst->plot = 1;

	inst->starting_time = second();

	inst->edge_weight = 1; // 0 ATT, 1 EUC_2D default

	inst->greedy_starts = 10;

	inst->vns_kmin = 1;
	inst->vns_kmax = 5;
	inst->vns_learning_rate = 1.0;
	inst->vns_jumps = 0;

	inst->tabu_min = 0;
	inst->tabu_max = 0;
	inst->tabu_tenure = 0;
	inst->tabu_noimprove = 0;

	inst->hard_fixing_percentage = 0.30;
	inst->hard_fixing_local_time = 15.0;

	inst->local_branch_k = 5;
	inst->local_branch_step = 5;

	strcpy(inst->csv_filename, "");
}

/**
 * Generate three distinct random cut indices i < j < k
 * for a 3-opt move on a tour.
 *
 * @param tsp Pointer to the TSP instance
 * @param i Pointer to store index i
 * @param j Pointer to store index j
 * @param k Pointer to store index k
 * @return void
 */
static void generate_three_opt_indices(instance *tsp, int *i, int *j, int *k)
{
	int nnodes = tsp->nnodes;
	while (1)
	{
		*i = rand() % (nnodes - 3);
		*j = *i + 1 + rand() % (nnodes - *i - 2);
		*k = *j + 1 + rand() % (nnodes - *j - 1);
		if (*k < nnodes - 1)
			break;
	}
}

/**
 * Copy a segment of the tour from src to dst. It works both in normal or reversed order, updating
 * the write position in the destination array.
 *
 * @param dst Destination tour array
 * @param src Source tour array
 * @param pos Pointer to current write position in dst
 * @param start Start index in src
 * @param end End index in src
 * @param reverse If 1, copy in reverse order; if 0, copy normally
 * @return void
 */
static void copy_segment(int *dst, int *src, int *pos, int start, int end, int reverse)
{
	if (reverse)
	{
		for (int i = end; i >= start; i--)
			dst[(*pos)++] = src[i];
	}
	else
	{
		for (int i = start; i <= end; i++)
			dst[(*pos)++] = src[i];
	}
}

/**
 * Apply a 3-opt move to the tour based on given type and cut indices.
 * The move reconnects the segments (i+1 to j), (j+1 to k), and (k+1 to end)
 * based on the specified 3-opt move type.
 *
 * @param tsp Pointer to the TSP instance
 * @param src_tour Source tour array
 * @param dst_tour Destination tour array to store the result
 * @param i First cut index
 * @param j Second cut index
 * @param k Third cut index
 * @param type Type of 3-opt move to apply (enum ThreeOptType)
 * @return void
 */
static void apply_three_opt_move(instance *tsp, int *src_tour, int *dst_tour, int i, int j, int k, enum ThreeOptType type)
{
	int pos = 0;
	int nnodes = tsp->nnodes;

	// Copy prefix 0..i
	memcpy(dst_tour, src_tour, sizeof(int) * (i + 1));
	pos = i + 1;

	int reverse_A = 0, reverse_B = 0, reverse_C = 0;

	switch (type)
	{
	case TYPE_0:
		break;
	case TYPE_1:
		reverse_A = 1;
		break;
	case TYPE_2:
		reverse_B = 1;
		break;
	case TYPE_3:
		reverse_C = 1;
		break;
	case TYPE_4:
		reverse_A = reverse_B = 1;
		break;
	case TYPE_5:
		reverse_A = reverse_C = 1;
		break;
	case TYPE_6:
		reverse_B = reverse_C = 1;
		break;
	}

	copy_segment(dst_tour, src_tour, &pos, i + 1, j, reverse_A);
	copy_segment(dst_tour, src_tour, &pos, j + 1, k, reverse_B);
	copy_segment(dst_tour, src_tour, &pos, k + 1, nnodes - 1, reverse_C);

	dst_tour[nnodes] = dst_tour[0]; // Close tour
}

/**
 * Apply the 3-opt heuristic to a solution.
 * Randomly selects a valid (i, j, k) triple and evaluates all 7 possible
 * 3-opt reconnection types. Applies the best one (lowest total cost).
 *
 * @param tsp Pointer to the TSP instance
 * @param sol Pointer to the solution to update
 * @return 0 if the 3-opt heuristic is applied successfully, -1 on error
 */
int apply_three_opt(instance *tsp, solution *sol)
{
	if (!tsp->cost_matrix || tsp->nnodes <= 0)
		return -1;

	int i, j, k;
	generate_three_opt_indices(tsp, &i, &j, &k);

	int nnodes = tsp->nnodes;
	int *new_tour = malloc(sizeof(int) * (nnodes + 1));
	int *best_tour = malloc(sizeof(int) * (nnodes + 1));
	if (!new_tour || !best_tour)
		return -1;

	double best_cost = CPX_INFBOUND;

	for (int move = 0; move < 7; move++)
	{
		apply_three_opt_move(tsp, sol->tour, new_tour, i, j, k, move);

		double cost = 0;
		for (int n = 0; n < nnodes; n++)
			cost += tsp->cost_matrix[new_tour[n]][new_tour[n + 1]];

		best_cost = cost;
		memcpy(best_tour, new_tour, sizeof(int) * (nnodes + 1));
	}

	memcpy(sol->tour, best_tour, sizeof(int) * (nnodes + 1));
	sol->cost = best_cost;

	free(new_tour);
	free(best_tour);

	return 0;
}

/**
 * Get the relative time weight of a single method, used to
 * partition the total time limit among multiple sub‐methods.
 *
 * @param method_name name of the method
 * @return the weight of the method
 */
double get_method_weight(const char *method_name)
{
	if (strcmp(method_name, "benders") == 0 || strcmp(method_name, "branch_and_cut") == 0 ||
		strcmp(method_name, "hard_fix") == 0 || strcmp(method_name, "local_branch") == 0)
		return 0.9;
	if (strcmp(method_name, "tabu") == 0)
		return 0.1;
	if (strcmp(method_name, "n_n") == 0 || strcmp(method_name, "random") == 0 || strcmp(method_name, "two_opt") == 0)
		return 0.1;
	if (strcmp(method_name, "vns") == 0)
		return 0.5;
	return 1.0; // default weight
}
