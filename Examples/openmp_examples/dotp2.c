#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N_MAX 100000
#define Max_thr 128 
double x[N_MAX],y[N_MAX];

int main() {
  int i, n = N_MAX;
  double t ,tx;
/*-------------------- set number of threads */
  omp_set_num_threads ( 4 ) ; //  set # threads
/*-------------------- divide dot product in equal parts*/
/*-------------------- generate random vectors */
  tx = 0.0;
  for (i = 0; i < n; i++){
    x[i] = (double) rand()/ ((double) RAND_MAX);
    y[i] = (double) rand()/ ((double) RAND_MAX);
/*-------------------- this computes exact ddot */
    tx += x[i]*y[i] ;
  }
  
  t = 0.0;
# pragma omp parallel for reduction (+:t)
  for ( i=0; i < n ; i ++) {
    t += x[i]*y[i];
  }
/*-------------------- print both results */
  printf("-- DotProduct  t = %f  ;   exact tx= %f\n",t,tx);
}
