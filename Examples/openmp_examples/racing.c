#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N_MAX 10000
int main() {
int i;
 double fx ,t, fsum =0.0;
#pragma omp parallel for
for (i = 1; i <= N_MAX; i++) {
  fx = (double)i ;
  fsum += fx;
 }
 t = (double) N_MAX*(N_MAX+1)/2.0;
 printf("-- sum %f   exact %f \n", fsum,t); 
}
