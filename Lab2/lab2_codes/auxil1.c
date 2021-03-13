#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define MAX_LINE 210

/*-------------------- begin prototyping */
int saxpy_(int *n, float *a, float *x, int *incx, float *y, int *incy);
int sscal_(int *n, float *alpha, float *x, int *inc);
void scopy_(int *nfeat, float *x, int *incx, float *y, int *incy);
void get_rand_ftr(float *ctr, float *fdata, int m, int nfeat);
/*-------------------- end prototyping */

void get_rand_ftr(float *ctr, float *fdata, int m, int nfeat){
// gets a random convex  combination of all samples
  float tot, t; 
  int j, one=1;
  /*-------------------- initialize to zero */
  for (j=0; j<nfeat; j++)
    ctr[j] = 0.0;
  tot = 0.0;
  /*-------------------- loop over all samples*/ 
  for (j=0; j<m; j++) {
    t = (float)(rand() / (float) RAND_MAX);
    t = t*t;
    //if (t < 0.5){
    //    for (k=0; k<nfeat; k++)      ctr[k] += t*fdata[j*nfeat+k];
    saxpy_(&nfeat,&t,&fdata[j*nfeat],&one,ctr,&one);
    tot +=t;    
  }
  tot = 1.0/tot;
  sscal_(&nfeat,&tot,ctr,&one);
}

int assign_ctrs(float *dist, int k){
  float min;
  min = dist[0];
  int i, ctr = 0;
  for(i=1; i<k; i++) {
    if(min > dist[i]) {
      ctr = i;
      min = dist[i];
    }
  }
  return ctr;
}

/*-------------------- reading data */
int read_csv_matrix(float *mtrx, char file_name[], int* nrow, int *nfeat){
/* -------------------- reads data from a csv file to mtrx */
  FILE *finputs; 
  char line[MAX_LINE], subline[MAX_LINE];
  
  if (NULL == (finputs = fopen(file_name, "r" )))
    exit(1);
  memset(line,0,MAX_LINE);
  //  
  int k, j, start, rlen, lrow=0, lfeat=0, jcol=0, jcol0=0, len, first= 1;
  char *delim;
  delim =",";
  /*-------------------- big while loop */
  while(fgets(line,MAX_LINE,finputs)) {
    if(first) {
//--------------------ignore first line of csv file 
      first = 0;
      continue;
    }
    len = strlen(line);
    lrow++;
    start = 0;
/*-------------------- go through the line */
    for (j=0; j<len; j++){
      if (line[j] == *delim || j==len-1) {
	k = j-start;
	memcpy(subline,&line[start],k*sizeof(char));
	//-------------------- select items to drop here --*/
	if (start >0){     //  SKIPPING THE FIRST RECORD 
	  subline[k] = '\0';
	  mtrx[jcol++] = atof(subline);
	}
	start = j+1;
      }
    }
/*-------------------- next row */
    rlen  = jcol -jcol0;
    jcol0 = jcol;
    if (lrow == 1) lfeat = rlen;
/*-------------------- inconsistent rows */
    if (rlen != lfeat)  return(1);
  }
/*-------------------- done */
  fclose(finputs);
  //  for (j=0; j<jcol; j++)    printf(" %e  \n",mtrx[j]);
  *nrow = lrow;
  *nfeat = lfeat;
  return(0);
}
/*-------------------- assign a center to an item */

float dist2(float *x, float *y, int len){
  int i;
  float dist = 0;
  for(i=0; i<len; i++){
    dist += (x[i] - y[i]) * (x[i] - y[i]);
  }
  return dist/len;
}

/*=======================================================================*/

int MyKmeans_p(float *fdata, int *clustId, int *counter, int *params,
	       float tol, MPI_Comm comm) {
/*==================================================
  IN: 
    fdata     = float* (nfeat*m)   = input data (local to this process).
    params[ ]  = int* contains 
    params[0] = Nc     = number of clusters
    params[1] = m      =  number of samples
    params[2] = nfeat  = number of features
    params[3] = maxits = max number of Kmeans iterations
    tol       = tolerance for determining if centers have converged.
    comm      = communicator
  OUT:
   clustId[i] = cluster of sample i for i=1...m
   counter[j] = size of cluster j for j=1:Nc
   ===================================================*/
  // some declarations
  int nprocs, myid;
  /*-------------------- START */
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&myid);
  //-------------------- unpack params.  
  int Nc = params[0];
  int m  = params[1]; 
  int nfeat = params[2]; 
  int maxits = params[3];
  //int NcNf =Nc*nfeat;
  /*-------------------- replace these by your function*/   
  float ctr[Nc][nfeat];
  int global_cnt, local_cnt;
  double global_sum, local_sum = 0.0;
  float dists[Nc];
  int i, j, k, total_cnt;
  total_cnt = 0;
  local_cnt = Nc;
  /*-------------------- Phase 0: initialization*/
  MPI_Allreduce(&local_cnt, &global_cnt, 1, MPI_INT, MPI_SUM, comm);
  //MPI_Barrier(comm);
  for (j=0; j<Nc; j++) {
    get_rand_ftr(ctr[j], fdata, m, nfeat); 
  }
  for (j=0; j<nfeat; j++) {
    for (i=0; i<Nc; i++) {
      local_sum = ctr[i][j];
      MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
      MPI_Barrier(comm);
      ctr[i][j] = (float) (global_sum / (double)global_cnt);
    }
  }
  //printf("Phase 0 complete.");
  //MPI_Barrier(comm);
  /*--------------------- Phase 1: get cluster ids*/
  while(total_cnt < maxits) {
  for (i=0; i<m; i++) {
    for (j=0; j<Nc; j++) {
      dists[j] = dist2(ctr[j], &fdata[i], nfeat);
    }
    clustId[i] = assign_ctrs(dists, Nc);
  }
  for (j=0; j<Nc; j++)
    counter[j] = 0.0;
  for (i=0; i<m; i++) {
    for (j=0; j<Nc; j++) {
      if (clustId[i] == j)
        counter[j]++;
    }
  }
  for (j=0; j<Nc; j++) {
    local_cnt = counter[j];
    MPI_Allreduce(&local_cnt, &global_cnt, 1, MPI_INT, MPI_SUM, comm);
    MPI_Barrier(comm);
    counter[j] = global_cnt;
  }
  //printf("Phase 1 complete.");
  //MPI_Barrier(comm);
  /*---------------------Phase 2: reset centoids*/
  for (j=0; j<Nc; j++) {
    //local_cnt = counter[j];
    //MPI_Allreduce(&local_cnt, &global_cnt, 1, MPI_INT, MPI_SUM, comm);
    for (i=0; i<nfeat; i++) {
      local_sum = 0.0;
      for (k=0; k<m; k++) {
	if (clustId[k] == j) {
          local_sum += (&fdata[k])[i];
	}
      }
      MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
      MPI_Barrier(comm);
      /*if (counter[j] == 0) {
        get_rand_ftr(ctr[j], fdata, m, nfeat);
	local_sum = ctr[j][i];
	MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
	local_cnt = Nc;
	MPI_Allreduce(&local_cnt, &global_cnt, 1, MPI_INT, MPI_SUM, comm);
      }*/
      if (counter[j] != 0)
        ctr[j][i] = (float) (global_sum / (double)counter[j]);
    }
  }
  //MPI_Barrier(comm);
  for (j=0; j<Nc; j++) {
    if (counter[j] == 0) {
      get_rand_ftr(ctr[j], fdata, m, nfeat);
      for (i=0; i<nfeat; i++) {
        local_sum = ctr[j][i];
	MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
	//MPI_Barrier(comm);
	local_cnt = Nc;
	MPI_Allreduce(&local_cnt, &global_cnt, 1, MPI_INT, MPI_SUM, comm);
	MPI_Barrier(comm);
	ctr[j][i] = (float) (global_sum / (double)global_cnt);
      }
    }
  }
  //MPI_Barrier(comm);
  //printf("Phase 2 complete.");
  total_cnt++;
  }
  return(0);
}

