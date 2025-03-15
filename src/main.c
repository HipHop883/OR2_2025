#include "../include/tsp.h"
#include "chrono.h"
//unset GTK_PATH to avoid conflict with gnuplot

int main(int argc, char **argv) 
{ 
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }       
	if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }

	double t1 = second(); 
	instance inst;
	
	parse_command_line(argc, argv, &inst);     
	
	//printf(" file %s has %d non-empty lines\n", inst.input_file, number_of_nonempty_lines(inst.input_file)); exit(1);

	read_input(&inst);  
	
	
	printf("Number of nodes: %d\n", inst.nnodes);
	//print_nodes(&inst);

	// Allocate memory for best_sol
    inst.best_sol = (double *) calloc(inst.nnodes + 1, sizeof(double));

	if(strcmp(inst.method, "n_n") == 0) {
		nearest_neighbor(&inst);			// Nearest neighbor heuristic
		check_time(&inst, t1);
	} else if(strcmp(inst.method, "n_n+two_opt") == 0) {
		nearest_neighbor(&inst);			// Nearest neighbor heuristic with two_opt
		two_opt(&inst);					
		check_time(&inst, t1);
	} else if(strcmp(inst.method, "random+two_opt") == 0) {
		random_path(inst.best_sol, inst.nnodes, inst.randomseed);	// Random path with two_opt
		two_opt(&inst);					
		check_time(&inst, t1);
	} else if(strcmp(inst.method, "random") == 0) {
		random_path(inst.best_sol, inst.nnodes, inst.randomseed);	// Random path
		check_time(&inst, t1);
	} else {
		print_error("USAGE: -method [nearest_neighbor|n_n+two_opt|random+two_opt|random]");
	}

	// Print the cost of the best solution
	cost_path(&inst);

	/* Print the best solution path
    printf("Best solution path:\n");
    print_path(&inst, inst.best_sol, inst.nnodes);
	*/

	double t2 = second(); 
    
	// Plot the solution path using Gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set xlabel 'X'\n");
        fprintf(gnuplotPipe, "set ylabel 'Y'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "set key top right\n"); // Enable legend and set position
        fprintf(gnuplotPipe, "plot '-' with linespoints lt rgb 'red' lw 2 pt 7 ps 1.5 title 'TSP-%s'\n", inst.method);
        for (int i = 0; i <= inst.nnodes; i++) {
            int node = (int)inst.best_sol[i];
            fprintf(gnuplotPipe, "%lf %lf\n", inst.xcoord[node], inst.ycoord[node]);
        }
        fprintf(gnuplotPipe, "e\n");
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        print_error("Error opening Gnuplot. Make sure Gnuplot is installed and in your PATH.");
    }

	if ( VERBOSE >= 1 )   
	{
		printf("... TSP problem solved in %lf sec.s\n", t2-t1);  
	}

	//free memory
	free_instance(&inst);
	return 0; 
}       
