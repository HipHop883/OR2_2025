#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <stdlib.h>

/**
 * Get current time in seconds
 * @return current time in seconds
 */
double myWallTime()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9 * ((double)ts.tv_nsec);
}

/**
 * Get current time in seconds
 * @return current time in seconds
 */
double second()
{
    double t = myWallTime();
    return t;
}
