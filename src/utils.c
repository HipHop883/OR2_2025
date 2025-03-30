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
    double t = (double)ts.tv_sec + 1.0e-9 * ((double)ts.tv_nsec);

    return t;
}

double rand01(int seed)
{
    srand(seed);
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
    return (*(int *)a - *(int *)b);
}
