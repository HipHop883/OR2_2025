#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <stdlib.h>
#include <string.h>

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
 * Check cmd line arguments for the --plot or -p flag
 * @param argc number of arguments
 * @param argv array of arguments
 * @return 1 if the --plot or -p flag is present, 0 otherwise
 */
int should_plot(int argc, char **argv)
{
    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--plot") || !strcmp(argv[i], "-p"))
        {
            return 1;
        }
    }
    return 0;
}
