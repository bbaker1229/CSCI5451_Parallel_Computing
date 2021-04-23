#include <mpi.h>
#include <math.h> 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "inc.h"

void OneStep (PointProb** myMatrix, double** Pij, double** old_Pij,
	      DomainPtr Dom,  MPI_Comm comm2D){
/*-------------------- Does one iteration */
  int i,j;
  MPI_Status status;
  int ni =Dom->ni, nj=Dom->nj;
  double *sendNS, *sendEW;
  sendNS = (double*) malloc(sizeof(double)*nj);
  sendEW = (double*) malloc(sizeof(double)*ni);

  // Initialize ghost vectors
  for(i=0;i<ni;i++) {
    Dom->west[i] = 0.0;
    Dom->east[i] = 0.0;
  }
  for(j=0;j<nj;j++) {
    Dom->north[j] = 0.0;
    Dom->south[j] = 0.0;
  }

  // Initialize new Pij matrix
  for(i=0; i<ni; i++)
    for(j=0; j<nj; j++)
      Pij[i][j] = 0.0;

  // Calculate Pij entries
  for (i=0; i<ni; i++) {
    for (j=0; j<nj; j++) {
      Pij[i][j]       += old_Pij[i][j] * myMatrix[i][j].prob[0];
      if(i == 0)
        Dom->north[j] += old_Pij[i][j] * myMatrix[i][j].prob[3];
      else
        Pij[i-1][j]   += old_Pij[i][j] * myMatrix[i][j].prob[3];
      if(j == 0)
        Dom->west[i]  += old_Pij[i][j] * myMatrix[i][j].prob[1];
      else
        Pij[i][j-1]   += old_Pij[i][j] * myMatrix[i][j].prob[1];
      if(i == ni-1)
        Dom->south[j] += old_Pij[i][j] * myMatrix[i][j].prob[4];
      else 
        Pij[i+1][j]   += old_Pij[i][j] * myMatrix[i][j].prob[4];
      if(j == nj-1)
        Dom->east[i]  += old_Pij[i][j] * myMatrix[i][j].prob[2];
      else
        Pij[i][j+1]   += old_Pij[i][j] * myMatrix[i][j].prob[2];
    }
  }

  //Send the ghost vectors back and add them to their new vectors
  if(Dom->qj % 2) { //if odd j rank - recv east vectors
    MPI_Send(Dom->west, ni, MPI_DOUBLE, Dom->westNB, 0, comm2D);
    MPI_Recv(sendEW, ni, MPI_DOUBLE, Dom->eastNB, 1, comm2D, &status);
    for(i=0;i<ni;i++)
      Pij[i][nj-1] += sendEW[i];
  } else {
    MPI_Recv(sendEW, ni, MPI_DOUBLE, Dom->eastNB, 0, comm2D, &status);
    for(i=0; i<ni; i++)
      Pij[i][nj-1] += sendEW[i];
    MPI_Send(Dom->west, ni, MPI_DOUBLE, Dom->westNB, 1, comm2D);
  }

  if(Dom->qj % 2) { //if odd j rank - recv west vectors
    MPI_Send(Dom->east, ni, MPI_DOUBLE, Dom->eastNB, 0, comm2D);
    MPI_Recv(sendEW, ni, MPI_DOUBLE, Dom->westNB, 1, comm2D, &status);
    for(i=0;i<ni;i++)
      Pij[i][0] += sendEW[i];
  } else {
    MPI_Recv(sendEW, ni, MPI_DOUBLE, Dom->westNB, 0, comm2D, &status);
    for(i=0; i<ni; i++)
      Pij[i][0] += sendEW[i];
    MPI_Send(Dom->east, ni, MPI_DOUBLE, Dom->eastNB, 1, comm2D);
  }

  if(Dom->qi % 2) { //if odd i rank - recv south vectors
    MPI_Send(Dom->north, nj, MPI_DOUBLE, Dom->northNB, 0, comm2D);
    MPI_Recv(sendNS, nj, MPI_DOUBLE, Dom->southNB, 1, comm2D, &status);
    for(j=0;j<nj;j++)
      Pij[ni-1][j] += sendNS[j];
  } else {
    MPI_Recv(sendNS, nj, MPI_DOUBLE, Dom->southNB, 0, comm2D, &status);
    for(j=0; j<nj; j++)
      Pij[ni-1][j] += sendNS[j];
    MPI_Send(Dom->north, nj, MPI_DOUBLE, Dom->northNB, 1, comm2D);
  }

  if(Dom->qi % 2) { //if odd i rank - recv north vectors
    MPI_Send(Dom->south, nj, MPI_DOUBLE, Dom->southNB, 0, comm2D);
    MPI_Recv(sendNS, nj, MPI_DOUBLE, Dom->northNB, 1, comm2D, &status);
    for(j=0;j<nj;j++)
      Pij[0][j] += sendNS[j];
  } else {
    MPI_Recv(sendNS, nj, MPI_DOUBLE, Dom->northNB, 0, comm2D, &status);
    for(j=0; j<nj; j++)
      Pij[0][j] += sendNS[j];
    MPI_Send(Dom->south, nj, MPI_DOUBLE, Dom->southNB, 1, comm2D);
  }
}


double ComputErr(double** x, double** y, int ni, int nj, MPI_Comm comm2D){
  double max=0;
  double temp;
  int i,j;
  /*local computation of max{abs(u_k+1(i)-u_k(i))}*/
  
  for(i=0;i<ni;i++){
    for(j=0;j<nj;j++){
      temp=fabs(x[i][j]-y[i][j]);
      if(temp>max)
        max=temp;
    }
  }
/*--------------------allreduce the results of max at root node*/
  MPI_Allreduce(&max,&temp,1,MPI_DOUBLE,MPI_MAX,comm2D);
  return temp;
}

int normliz(double** x, int ni, int nj, MPI_Comm comm2D){
  double temp=0.0;
  int i,j;
  /*local computation of max{abs(u_k+1(i)-u_k(i))}*/
  for(i=0;i<ni;i++)
    for(j=0;j<nj;j++)
      temp+=x[i][j];
/*--------------------allreduce the results of max at root node*/
  MPI_Allreduce(&temp,&temp,1,MPI_DOUBLE,MPI_SUM,comm2D);
  if (temp == 0.0)
    return 1;
  for(i=0;i<ni;i++)
    for(j=0;j<nj;j++)
      x[i][j] /=temp;   
  return 0;
}


void get_prob(int i, int j, int ni, int nj, int pi, int pj,
		 int *dims, double *prob) {
/* get probability stencil */
/* ------------------------*/
/* i, j = local coordinates of point in each cell 
   ni, nj = dimensions of cell [assumed to be all the same size]
   pi, pj = process id in cartesian topology 
   dims   = array of size to containing the max number of 
            processes in each direction of cartesian mesh 
------------------------------------------------------------*/
  double tw, te, tn, ts;
  int ii, jj, niT, njT;
  /*-------------------- global indices */
  ii = pi*ni + i;
  jj = pj*nj + j;
  /*-------------------- global matrix dimensions */
  niT = dims[0]*ni;
  njT = dims[1]*nj;  
  // T1: pl = .15  pr = .3   == all multiplied by .25 times fun. of i,j
  // T2  p = p = .2
  te = 0.25* 0.30 * (double)(njT - jj+1)/ (double) njT;
  tw = 0.25* 0.15 * (double)(njT + jj-1)/ (double) njT;
  ts = 0.25* 0.20 * (double)(niT - ii+1)/ (double) niT;
  tn = 0.25* 0.20 * (double)(niT + ii-1)/ (double) niT;
  /*-------------------- save trans. probabilities */
  prob[1] = tw;      // west
  prob[2] = te;      // east
  prob[3] = tn;      // north 
  prob[4] = ts;      // south
  prob[0] = 1.0 - (te+tw+tn+ts);
}
