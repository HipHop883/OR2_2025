#include "tsp.h"
#include "tsp_greedy.h"
#include "utils.h"

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		printf("Usage: %s -help for help\n", argv[0]);
		return EXIT_FAILURE;
	}

	if (VERBOSE >= 2)
	{
		for (int a = 0; a < argc; a++)
			printf("%s ", argv[a]);
		printf("\n");
	}

	int max_runs = runs(argc, argv);

	double *best_costs = (double *)malloc(max_runs * sizeof(double));
	if (!best_costs)
	{
		print_error("Memory allocation failed for best_costs array");
		return EXIT_FAILURE;
	}

	double *execution_times = (double *)malloc(max_runs * sizeof(double));
	if (!execution_times)
	{
		print_error("Memory allocation failed for execution_times array");
		return EXIT_FAILURE;
	}

	instance inst;
	init(&inst);
	solution sol;

	if (parse_command_line(argc, argv, &inst))
	{
		print_error("Error parsing command line");
		return EXIT_FAILURE;
	}

	int master_seed = inst.randomseed;

	for (int i = 0; i < max_runs; i++)
	{
		set_seed(master_seed + i);

		if (load_instance(&inst))
		{
			print_error("Error reading input");
			free_instance(&inst);
			free_sol(&sol);

			return EXIT_FAILURE;
		}

		printf("Number of nodes: %d\n", inst.nnodes);

		double t1 = second();

		if (execute_selected_method(&inst, &sol))
		{
			print_error("Error running method\n");
			free_instance(&inst);
			free_sol(&sol);

			return EXIT_FAILURE;
		}

		double elapsed_time = second() - t1;

		printf("Total tour cost: %.2lf\n", sol.cost);

		execution_times[i] = elapsed_time;
		best_costs[i] = sol.cost;

		// Plot the solution path using Gnuplot
		if (inst.plot == 0)
		{
			FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
			if (gnuplotPipe)
			{
				char filename[256];
				sprintf(filename, "plot/TSP_%s.png", inst.method);
				fprintf(gnuplotPipe, "set terminal png\n");
				fprintf(gnuplotPipe, "set output '%s'\n", filename);
				fprintf(gnuplotPipe, "set xlabel 'X'\n");
				fprintf(gnuplotPipe, "set ylabel 'Y'\n");
				fprintf(gnuplotPipe, "set grid\n");
				fprintf(gnuplotPipe, "set key top right\n");
				fprintf(gnuplotPipe, "set termoption noenhanced\n");
				fprintf(gnuplotPipe, "plot '-' with linespoints lt rgb 'red' lw 2 pt 7 ps 1.5 title 'TSP-%s (Cost: %.2lf)'\n", inst.method, sol.cost);
				for (int i = 0; i <= inst.nnodes; i++)
				{
					int node = sol.tour[i];
					fprintf(gnuplotPipe, "%lf %lf\n", inst.xcoord[node], inst.ycoord[node]);
				}
				fprintf(gnuplotPipe, "e\n");
				fflush(gnuplotPipe);
				pclose(gnuplotPipe);
				printf("Plot saved in: %s\n", filename);
			}
			else
			{
				print_error("Error opening Gnuplot. Make sure Gnuplot is installed and in your PATH.");
			}
		}

		if (VERBOSE >= 1)
		{
			printf("TSP solved in %lf seconds\n", elapsed_time);
		}

		free_sol(&sol);
	}

	free_instance(&inst);

	if (!strcmp(inst.method, "benders") || !strcmp(inst.method, "branch_and_cut") ||
		!strcmp(inst.method, "tabu+benders") || !strcmp(inst.method, "tabu+branch_and_cut") ||
		!strcmp(inst.method, "two_opt+benders") || !strcmp(inst.method, "two_opt+branch_and_cut"))
	{
		update_perf_csv(&inst, execution_times, max_runs);
	}
	else
	{
		update_perf_csv(&inst, best_costs, max_runs);
	}

	free(best_costs);
	free(execution_times);

	return EXIT_SUCCESS;
}
