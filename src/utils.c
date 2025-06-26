#define _POSIX_C_SOURCE 200809L
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "utils.h"

/**
 * Get current time in seconds
 * @return current time in seconds
 */
double second()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9 * ts.tv_nsec;
}
/**
 * Set the random seed and run rome random numbers to improve randomness
 * @param seed the seed to set
 */
void set_seed(int seed)
{
    srand(seed);
    // To imporve randomness
    for (int i = 0; i < 10000; i++)
        rand();
}

/**
 * Generate a random number between 0 and 1
 * @return a random number between 0 and 1
 */
double rand01()
{
    return (double)rand() / RAND_MAX;
}

/**
 * Compare two integers for qsort
 * @param a integer a
 * @param b integer b
 * @return the difference between a and b
 */
int compar(const void *a, const void *b)
{
    int av = *(int *)a;
    int bv = *(int *)b;
    return (av > bv) - (av < bv);
}
/**
 * Get the number of runs to perform from the command line arguments
 * @param argc number of arguments
 * @param argv array of arguments
 * @return the number of runs, default is 1
 */
int runs(int argc, char **argv)
{
    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--runs") || !strcmp(argv[i], "-rs"))
        {
            if (i + 1 < argc)
            {
                return atoi(argv[++i]);
            }
        }
    }

    return 1;
}

/**
 * Update the performance CSV file for the current configuration.
 *
 * The CSV file is expected to have:
 *  - A header line: <num_alg>,<instance_id1>,<instance_id2>,...,<instance_idN>
 *  - Each subsequent line: runX,<cost1>,<cost2>,...,<costN>
 *
 * If the current instance_id (constructed from inst's VNS parameters)
 * is not present in the header, it is added (and the number of algorithms is increased).
 * Then, for each run (from 1 to num_runs), the cost value in the corresponding column
 * is updated with run_results[i]. If a run row doesn't exist, it is appended.
 *
 * @param inst       The current instance (provides method, timelimit, VNS params, seed, etc.)
 * @param run_results  An array of double values, one per run (length num_runs).
 * @param num_runs   Number of runs in the current execution.
 */
void update_perf_csv(const instance *inst, double *run_results, int num_runs)
{
    char csv_filename[256];
    sprintf(csv_filename, "runs/%s_%.2lf.csv", inst->method, inst->timelimit);

    // Build the instance identifier (column header) based on VNS parameters.
    char instance_id[256];
    if (!strcmp(inst->method, "vns"))
    {
        if (inst->vns_jumps > 0)
        {
            sprintf(instance_id, "fixed_kick_%d_lr_%.2lf_seed_%d",
                    inst->vns_jumps, inst->vns_learning_rate, inst->randomseed);
        }
        else
        {
            sprintf(instance_id, "rand_kick_%d_%d_lr_%.2lf_seed_%d",
                    inst->vns_kmin, inst->vns_kmax, inst->vns_learning_rate, inst->randomseed);
        }
    }
    else if (!strcmp(inst->method, "tabu"))
    {
        sprintf(instance_id, "tabu_tenure_%d_%d_%d_noimpr_%d_seed_%d",
                inst->tabu_min, inst->tabu_tenure, inst->tabu_max,
                inst->tabu_noimprove, inst->randomseed);
    }
    else if (!strcmp(inst->method, "benders"))
    {
        sprintf(instance_id, "benders_seed_%d", inst->randomseed);
    }
    else if (!strcmp(inst->method, "branch_and_cut"))
    {
        sprintf(instance_id, "branch_and_cut_seed_%d", inst->randomseed);
    }
    else
    {
        // TODO other methods
        strcpy(instance_id, inst->method);
    }

    FILE *fp = fopen(csv_filename, "r");
    char **lines = NULL;
    int line_count = 0;
    if (fp != NULL)
    {
        char buffer[1024];
        while (fgets(buffer, sizeof(buffer), fp))
        {
            buffer[strcspn(buffer, "\n")] = '\0'; // remove newline
            lines = realloc(lines, sizeof(char *) * (line_count + 1));
            lines[line_count] = strdup(buffer);
            line_count++;
        }
        fclose(fp);
    }

    // Process header.
    char header[1024] = "";
    int num_cols = 0;   // current number of algorithms (columns)
    int col_index = -1; // column index for instance_id (1-indexed, as header fields 1..n)
    if (line_count > 0)
    {
        strcpy(header, lines[0]);

        char header_copy[1024];
        strcpy(header_copy, header);
        char *token;
        char *saveptr;
        token = strtok_r(header_copy, ",", &saveptr);
        if (token)
            num_cols = atoi(token);
        else
            num_cols = 0;

        int index = 1;
        char *found = NULL;
        token = strtok_r(NULL, ",", &saveptr);
        while (token != NULL)
        {
            if (strcmp(token, instance_id) == 0)
            {
                found = token;
                col_index = index;
                break;
            }
            index++;
            token = strtok_r(NULL, ",", &saveptr);
        }
        if (!found)
        {
            num_cols++;
            col_index = num_cols;
            char new_header[1024] = "";
            sprintf(new_header, "%d", num_cols);

            char *ptr = strchr(header, ',');
            if (ptr != NULL)
            {
                strcat(new_header, ptr);
            }

            strcat(new_header, ",");
            strcat(new_header, instance_id);
            free(lines[0]);
            lines[0] = strdup(new_header);
        }
    }
    else
    {

        num_cols = 1;
        col_index = 1;
        char new_header[1024] = "";
        sprintf(new_header, "1,%s", instance_id);
        lines = realloc(lines, sizeof(char *) * (line_count + 1));
        lines[line_count] = strdup(new_header);
        line_count++;
    }

    // Each run row should have (num_cols + 1) fields: first is the run label, then one cost per algorithm.
    // For each run i (i from 0 to num_runs - 1), update the cost at column col_index.
    for (int i = 0; i < num_runs; i++)
    {
        int row_index = i + 1; // header is row 0
        char new_row[1024] = "";
        if (row_index < line_count)
        {
            char row_copy[1024];
            strcpy(row_copy, lines[row_index]);

            char *tokens[128];
            int token_count = 0;
            char *tok;
            char *saveptr2;
            tok = strtok_r(row_copy, ",", &saveptr2);
            while (tok != NULL && token_count < 128)
            {
                tokens[token_count++] = tok;
                tok = strtok_r(NULL, ",", &saveptr2);
            }
            if (token_count > 0) // Check if at least one token was found
            {
                sprintf(new_row, "%s", tokens[0]);
                for (int j = 1; j <= num_cols; j++)
                {
                    strcat(new_row, ",");
                    if (j == col_index)
                    {
                        char cost_str[64];
                        sprintf(cost_str, "%.6lf", run_results[i]);
                        strcat(new_row, cost_str);
                    }
                    else if (j < token_count)
                    {
                        strcat(new_row, tokens[j]);
                    }
                    else
                    {
                        strcat(new_row, "");
                    }
                }
                free(lines[row_index]);
                lines[row_index] = strdup(new_row);
            }
            else
            {
                // Handle the case where the row is empty or has no tokens.
                // Initialize new_row differently in this case.
                sprintf(new_row, "run%d", i + 1);
                for (int j = 1; j <= num_cols; j++)
                {
                    strcat(new_row, ",");
                    if (j == col_index)
                    {
                        char cost_str[64];
                        sprintf(cost_str, "%.2lf", run_results[i]);
                        strcat(new_row, cost_str);
                    }
                    else
                    {
                        strcat(new_row, "");
                    }
                }
                free(lines[row_index]);
                lines[row_index] = strdup(new_row);
            }
        }
        else
        {
            sprintf(new_row, "run%d", i + 1);
            for (int j = 1; j <= num_cols; j++)
            {
                strcat(new_row, ",");
                if (j == col_index)
                {
                    char cost_str[64];
                    sprintf(cost_str, "%.2lf", run_results[i]);
                    strcat(new_row, cost_str);
                }
                else
                {
                    strcat(new_row, "");
                }
            }
            lines = realloc(lines, sizeof(char *) * (line_count + 1));
            lines[line_count] = strdup(new_row);
            line_count++;
        }
    }

    fp = fopen(csv_filename, "w");
    if (fp == NULL)
    {
        print_error("Failed to open CSV file for writing.");
        for (int i = 0; i < line_count; i++)
            free(lines[i]);
        free(lines);
        return;
    }
    for (int i = 0; i < line_count; i++)
    {
        fprintf(fp, "%s\n", lines[i]);
        free(lines[i]);
    }
    free(lines);
    fclose(fp);
}
