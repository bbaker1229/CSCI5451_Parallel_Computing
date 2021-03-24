#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#define MAX_LENGTH 100
/*-----------------------------------------------------------------------
 * this is to compute the exponential of a k numbers in a pipelined 
 * fashion -- set of values in this case. 
 * ----------------------------------------------------------------------*/
int main(int argc, char *argv[]){
  int myPE, East, West, n;
  double x[MAX_LENGTH], val;
  int k=10, j; 
  MPI_Status stat;
  MPI_Request req;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&n);
  MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
/*------------------- generate k numbers in [1 2]  */
  if (myPE == 0){
    for (j=0; j<k; j++)
      x[j] = 1.0 + ((double) j / (double) (k-1));
  }
  /*-------------------- broadcast value x from node 0*/
  MPI_Bcast(x, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
  East = myPE+1; West = myPE-1;
  /*-------------------- start */
  for (j=0; j<k; j++){
    if (myPE == n-1) {
      val = 1.0 + x[j] /( (double) myPE+1 ) ;
      MPI_Send(&val, 1, MPI_DOUBLE, West, j, MPI_COMM_WORLD) ;
    }   else    {
      MPI_Recv (&val, 1, MPI_DOUBLE, East, j, MPI_COMM_WORLD,&stat); 
      val = 1.0 + val *( x[j] /( (double) myPE+1));
      if (myPE>0) 
	MPI_Send(&val, 1, MPI_DOUBLE, West, j, MPI_COMM_WORLD) ;
      else
	printf(" x[j] = %e   val = %e \n", x[j],val);
    }
  }
  MPI_Finalize(); 
  return(0);
}

