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
  int len;
  int i, j, k, val, total_cnt = 0;
  int *global_cnt;
  float dists[Nc];
  float *ctr, *global_ctr, *local_sum, *oneline, *global_one;

  len = Nc*nfeat;
  ctr = (float *) malloc(len*sizeof(float));
  global_ctr = (float *) malloc(len*sizeof(float));
  local_sum = (float *) malloc(len*sizeof(float));
  oneline = (float *) malloc(nfeat*sizeof(float));
  global_one = (float *) malloc(nfeat*sizeof(float));
  global_cnt = (int *) malloc(Nc*sizeof(int));
  /*-------------------- Phase 0: initialization*/
  // Get random center values for each node
  for (k=0; k<Nc; k++)
    get_rand_ftr(&ctr[k*nfeat], fdata, m, nfeat); 
  // Find averages of all centers and sync to all nodes
  MPI_Allreduce(ctr, global_ctr, len, MPI_FLOAT, MPI_SUM, comm);
  for (k=0; k<Nc; k++)
    for (j=0; j<nfeat; j++)
      ctr[k*nfeat+j] = global_ctr[k*nfeat+j]/((float)nprocs); 
  /*--------------------- Phase 1: get cluster ids*/
  while(total_cnt < maxits) {
  // for each line calculate the distances to each center value
  // assign k value to line based on min distance.
  for (k=0; k<Nc; k++) {
    counter[k] = 0;
    for (j=0; j<nfeat; j++)
      local_sum[k*nfeat+j] = 0.0;
  }
  for (i=0; i<m; i++) {
    for (k=0; k<Nc; k++)
      dists[k] = dist2(&ctr[k*nfeat], &fdata[i*nfeat], nfeat);
    val = assign_ctrs(dists, Nc);
    clustId[i] = val;
    counter[val]++;
    for(j=0; j<nfeat; j++)
      local_sum[val*nfeat+j] += fdata[i*nfeat+j];
  }
  MPI_Allreduce(counter, global_cnt, Nc, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(local_sum, global_ctr, len, MPI_FLOAT, MPI_SUM, comm);
  for(k=0; k<Nc; k++)
    counter[k] = global_cnt[k]; 
  /*---------------------Phase 2: reset centoids*/
  // for each k value find the new center mean for each feature
  // and sync to all nodes if the counter value is not zero.
  for (k=0; k<Nc; k++) {
    if (counter[k] != 0) {
      for(j=0; j<nfeat; j++)
        ctr[k*nfeat+j] = global_ctr[k*nfeat+j]/((float)counter[k]);
    } 
    else {
      get_rand_ftr(&ctr[k*nfeat], fdata, m, nfeat);
      for(j=0; j<nfeat; j++)
        oneline[j] = ctr[k*nfeat+j];
      MPI_Allreduce(oneline, global_one, nfeat, MPI_FLOAT, MPI_SUM, comm);
      for(j=0; j<nfeat; j++)
        ctr[k*nfeat+j] = global_one[j]/((float)nprocs);
    }
  }
  // Increment the loop counter
  total_cnt++;
  }
  return(0);
}

