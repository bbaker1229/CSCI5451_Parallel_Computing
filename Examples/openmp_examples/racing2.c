#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N_MAX 100
int main() {
int i;
double fsum =0.0;
#pragma omp parallel for reduction (+:fsum)
for (i = 1; i <= N_MAX; i++) {
  fsum += (double) i;
 }
 printf("-- sum %f \n", fsum);
}
