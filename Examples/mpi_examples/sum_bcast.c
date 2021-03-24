#include "mpi.h"
#include <stdio.h>
#include <stdlib.h> 
#define MAXSIZE 1000

int main(int argc, char *argv[]){
  int n, myid, studt, numprocs;
  int data[MAXSIZE], locsum, sum, sum0; 
  int j, m, low, high; 

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
/*------------------------------------------------------------ */
/*-------------------- generate n random numbers ------------- */
  n = 100 ;
  sum0 = 0; 
/*---------- PUT THE a 3- DIGIT number here */ 
  studt = 231;   

  if (myid == 0) {
    fprintf(stdout," -- running partial sums -- output from proc 0\n");
    for (j=0; j<n; j++) {
/*-------------------- generate data and exact sum*/      
      data[j] = ((j+1)*117 + studt)  %  n;
      sum0 += data[j] ;
    }
  } 
/*------------------------------------------------------------ */
/*-------------------- broadcast data to processors ---------- */
  //MPI_Bcast(buff, len, datType, root,comm);
    MPI_Bcast(data, n, MPI_INT, 0, MPI_COMM_WORLD);
/*-------------------- calculate subsums  -------------------- */
  m = (int) ((n+numprocs-1) / numprocs); 
  low = myid * m; 
  high = low + m;
  if (high > n) high=n;
  locsum = 0;
  for (j=low; j<high; j++)
    locsum += data[j] ;
  /* printf("I got a sum of %d in %d\n", locsum, myid); */ 
  /* printf(" doing sum %d to %d\n", low, high); */
/*------------------------------------------------------------ */
/*-------------------- calculate sum of subsums (reduction)--- */
  //MPI_Reduce(sendbuf,recvbuff, cnt, type, OPE, root, comm);
  MPI_Allreduce(MPI_IN_PLACE, &locsum,1,MPI_INT,MPI_SUM, MPI_COMM_WORLD);
  if (myid == 0){ 
    printf("  ***** The total parallel computed sum is  %d \n", locsum);
    printf("  ***** The actual sum is                   %d \n", sum0);
  }
  MPI_Finalize(); 
}


