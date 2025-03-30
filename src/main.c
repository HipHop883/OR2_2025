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

	for (int i = 0; i < max_runs; i++)
	{
		instance inst;
		init(&inst);
		solution sol; // Safe init

		double t1 = second();

		// Parse command line
		if (parse_command_line(argc, argv, &inst))
		{
			print_error("Error parsing command line");
			return EXIT_FAILURE;
		}

		printf("Loading from file..\n\n\n\n");
		// Load problem instance
		if (load_instance(&inst))
		{
			printf("ERRORROROR");

			print_error("Error reading input");
			free_instance(&inst, &sol);
			return EXIT_FAILURE;
		}

		printf("Number of nodes: %d\n", inst.nnodes);

		// Allocate memory for the solution tour
		sol.tour = (int *)calloc(inst.nnodes + 1, sizeof(int));
		if (!sol.tour)
		{
			print_error("Memory allocation failed for solution");
			free_instance(&inst, &sol);
			return EXIT_FAILURE;
		}

		// Solve using the selected method
		if (execute_selected_method(&inst, &sol))
		{
			print_error("Error running method");
			free_instance(&inst, &sol);
			return EXIT_FAILURE;
		}

		// Print final cost
		printf("Total tour cost: %.2lf\n", sol.cost);

		double t2 = second();

		// Optional: plot the solution if --plot or -p is set
		if (should_plot(argc, argv))
		{
			FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
			if (gnuplotPipe)
			{
				fprintf(gnuplotPipe, "set xlabel 'X'\n");
				fprintf(gnuplotPipe, "set ylabel 'Y'\n");
				fprintf(gnuplotPipe, "set grid\n");
				fprintf(gnuplotPipe, "set key top right\n");
				fprintf(gnuplotPipe, "plot '-' with linespoints lt rgb 'red' lw 2 pt 7 ps 1.5 title 'TSP-%s'\n", inst.method);
				for (int i = 0; i <= inst.nnodes; i++)
				{
					int node = sol.tour[i];
					fprintf(gnuplotPipe, "%lf %lf\n", inst.xcoord[node], inst.ycoord[node]);
				}
				fprintf(gnuplotPipe, "e\n");
				fflush(gnuplotPipe);
				pclose(gnuplotPipe);
			}
			else
			{
				print_error("Gnuplot error. Make sure it's installed and in your PATH.");
			}
		}

		if (VERBOSE >= 1)
		{
			printf("TSP solved in %.2lf seconds\n", t2 - inst.starting_time);
		}

		free_instance(&inst, &sol);
	}

	return EXIT_SUCCESS;
}
