#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <unistd.h>
#define N_MAX 50000000 /* 50,000,000 */
int main(int argc, char* argv[])
{
#pragma omp parallel  
  printf("Hello, world.\n");
  return 0;
}
