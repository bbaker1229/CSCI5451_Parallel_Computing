#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N_MAX 100

int main() {
  int i, n=N_MAX;
  double t, fmax =0.0;
#pragma omp parallel for reduction (max:fmax) private (t)
  for (i = 1; i <= n; i++) {
    t = (double) n*i-i*i;
    if (fmax < t) 
      fmax = t; 
  }
  printf("-- max %f \n", fmax);
}

