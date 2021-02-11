#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N_MAX 100000
#define Max_thr 128 
double x[N_MAX],y[N_MAX];
double tt[Max_thr];

int main() {
  int i, i1, i2, it, m, nt, n = N_MAX;
  double t ,tx;
  nt = 4;
/*-------------------- set number of threads */
  omp_set_num_threads ( nt ) ; // nt = # threads
/*-------------------- divide dot product in equal parts*/
  m = (int) (1 + (n-1) / nt) ;      
/*-------------------- generate random vectors */
  tx = 0.0;
  for (i = 0; i < n; i++){
    x[i] = (double) rand()/ ((double) RAND_MAX);
    y[i] = (double) rand()/ ((double) RAND_MAX);
/*-------------------- this computes exact ddot */
    tx += x[i]*y[i] ;
  }
# pragma omp parallel for private (t, i1, i2, i)
  for ( it = 0; it < nt ; it ++) {
    i1 = it * m ;
    i2 = i1 + m ;
    if ( i2 > n ) i2 = n ;
    t = 0.0;
    for ( i = i1 ; i < i2 ; i ++ )
      t += x [ i ]* y [ i ];
    tt [it] = t ;
  }
/*-------------------- now sum partial dot products */ 
  t = 0.0;
  for (it = 0; it < nt ; it ++)
    t += tt [ it ];
/*-------------------- print both results */
  printf("-- DotProduct  t = %f  ;   exact tx= %f\n",t,tx);
}
