#include <sys/time.h>
#include <ctime>

double uni10_wtime()
{
  struct timeval time;
  return gettimeofday(&time, NULL) ? 0
                                   : (double)time.tv_sec + (double)time.tv_usec * 1e-6;
}

double uni10_cputime()
{
  return (double)clock() / CLOCKS_PER_SEC;
}
