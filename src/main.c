#include "tsp.h"
#include "tsp_greedy.h"
#include "utils.h"

// unset GTK_PATH to avoid conflict with gnuplot

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		printf("Usage: %s -help for help\n", argv[0]);
		exit(1);
	}

	if (VERBOSE >= 2)
	{
		for (int a = 0; a < argc; a++)
			printf("%s ", argv[a]);
		printf("\n");
	}

	double t1 = second();
	instance inst;
	init(&inst);
	solution sol;

	if (parse_command_line(argc, argv, &inst))
	{
		print_error("Error parsing command line");
		return -1;
	}

	if (load_instance(&inst))
	{
		print_error("Error reading input");
		free_instance(&inst, &sol);
		return -1;
	}

	printf("Number of nodes: %d\n", inst.nnodes);

	// Run the method and print the cost of the solution
	if (run_method(&inst, &sol))
	{
		print_error("Error running method\n");
		free_instance(&inst, &sol);
		return -1;
	}

	printf("The total cost is: %lf\n", sol.cost);

	/* Print the best solution path
	printf("Best solution path:\n");
	print_path(&inst, inst.best_sol, inst.nnodes);
	*/

	double t2 = second();

	// Plot the solution path using Gnuplot
	if (inst.plot == 0)
	{
		FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
		if (gnuplotPipe)
		{
			char filename[256];
			sprintf(filename, "plot/TSP_%s.png", inst.method); 	// Create the file name (es. TSP_2opt.png)
	
			fprintf(gnuplotPipe, "set terminal png\n");
			fprintf(gnuplotPipe, "set output '%s'\n", filename);
			fprintf(gnuplotPipe, "set xlabel 'X'\n");
			fprintf(gnuplotPipe, "set ylabel 'Y'\n");
			fprintf(gnuplotPipe, "set grid\n");
			fprintf(gnuplotPipe, "set key top right\n"); 		// Enable legend and set position
			fprintf(gnuplotPipe, "plot '-' with linespoints lt rgb 'red' lw 2 pt 7 ps 1.5 title 'TSP-%s'\n", inst.method);
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
		printf("... TSP problem solved in %lf sec.s\n", t2 - inst.starting_time);
	}

	// free memory
	free_instance(&inst, &sol);
	return 0;
}
