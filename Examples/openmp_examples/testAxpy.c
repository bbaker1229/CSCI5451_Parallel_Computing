#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h>

/*                   define MAX 50000000   */
#define N_MAX 1000000      /* 10^5 */
#define NITER 100       /* 100 iterations */

double x[N_MAX],y[N_MAX];

void* main ( )
{
  double t1, t2, td, target_work;
  int i, iter, nt, n;
  int j, k;  
  n = N_MAX;
  double wctime() ;
  float nops, alph;
/*-------------------- if niter entered online read it */

  printf(" Number of processors   %d\n",omp_get_num_procs()) ;
  printf(" Max number of threads  %d\n",omp_get_max_threads()); 
/*-------------------- loop over number of threads */
  for (nt = 1; nt <=128; nt*=2){
/*-------------------- initialize vectors to constants          */
  for (i = 0; i<n; i++){
    x[i] = 1.5;
    y[i] = 3.5;
  } 
/*-------------------- set number of threads */
  alph = 1.0/ (float)(NITER);
/*-------------------- * time iterations    */
  t1 = wctime();
  for (iter = 0; iter< NITER; iter++){
    omp_set_num_threads(nt);  
    /*#pragma omp parallel for schedule(guided,100)  */
    /*#pragma omp parallel for schedule(dynamic,10000)*/
#pragma omp parallel for schedule(runtime)
    for (i = 0; i < n-1; i++) {
      y[i] += alph*x[i];
    }
  }
  t2 = wctime();
  td = (t2 - t1)*1.e+06;
  nops = (float) NITER*n*2;
  printf(" ** num_threads = %4d, Mflops = %8.2f\n",nt,nops/td );
 }

}

