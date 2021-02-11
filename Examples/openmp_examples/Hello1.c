#include <stdio.h>
int omp_get_thread_num();
int omp_set_num_threads();
int omp_get_num_threads();
int main(){
  int i=12; int j;
  /* omp_set_num_threads(i); */
  
  omp_set_num_threads(i);
# pragma omp parallel 
  {
    j = omp_get_num_threads();
    printf("Thread number: %d -- %d \n", omp_get_thread_num(),j);
  }
  printf(" <<-- Out of threads \n"); 
}
