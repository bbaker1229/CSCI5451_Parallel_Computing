#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define BSIZE  16
#define NN 4000;

__global__ void MatAdd(int N, float *A, float *B, float *C){ 
   int i = blockIdx.x * blockDim.x + threadIdx.x;
   int j = blockIdx.y * blockDim.y + threadIdx.y;
   if (i < N && j < N)
       C[i*N+j] = A[i*N+j] + B[i*N+j]; 
} 

 void err_exit(char *message);
 float mat_add_check(int n, float *x, float *y, float *z)  {
 float s=0.0, t = 0.0, td = 0.0;
 for (int i=0; i<n; i++) {
       s  = y[i]+x[i]-z[i]; 
       t += s*s ;
       td += (x[i]*x[i]+y[i]*y[i]);
 }    
//-------------------- matrices are both zero
 if (td == 0.0) return(-1);
    else
//-------------------- normal return
   return(sqrt(s/td));
} 

int main() {
float *Ad, *Bd, *Cd; 
float  *A,  *B,  *C; 
int N, i, j; 
size_t MatSize;
float s;
//-------------------- set dimension N
 N = NN;

 char LineG[] = "Error allocating GPU  memory";
 char LineH[] = "Error allocating Host memory";

  
 MatSize = N*N*sizeof(float);
//-------------------- allocate on cpu
 A = (float *)malloc(MatSize);        
 B = (float *)malloc(MatSize);        
 C = (float *)malloc(MatSize);    
 if ((A==NULL) | (B==NULL) | (C==NULL) ) 
          err_exit(LineH);
//-------------------- allocate on GPU
 if (cudaMalloc((void **) &Ad, MatSize) != cudaSuccess) 
       err_exit(LineG);
 if (cudaMalloc((void **) &Bd, MatSize) != cudaSuccess) 
       err_exit(LineG);
 if (cudaMalloc((void **) &Cd, MatSize) != cudaSuccess) 
       err_exit(LineG);
//-------------------- fill arrays A,B

 for (i=0; i<N; i++) 
    for (j=0; j<N; j++) {
      A[i*N+j] = (float) rand() / (float) rand();
      B[i*N+j] = (float) rand() / (float) rand();
} 
//
//-------------------- copy matrices A,B+ to GPU memory
cudaMemcpy(Ad, A, MatSize, cudaMemcpyHostToDevice);
cudaMemcpy(Bd, B, MatSize, cudaMemcpyHostToDevice);
//-------------------- Kernel invocation
   dim3 dimBlock(BSIZE, 256/BSIZE);
   dim3 dimGrid((N + dimBlock.x-1) / dimBlock.x,
                (N + dimBlock.y-1) / dimBlock.y);
   MatAdd<<<dimGrid, dimBlock>>>(N, Ad, Bd, Cd);
//-------------------- see if things did execute 
 cudaError_t error = cudaGetLastError();
 if (error) {
     printf("CUDA error: %s \n",cudaGetErrorString(error));
     exit(1);
 }
//-------------------- Transfer result from GPU to CPU
cudaMemcpy(C, Cd, MatSize, cudaMemcpyDeviceToHost);
//-------------------- check whether addition was correct
s =  mat_add_check(N*N,A,B,C);
 
printf(" Mat dim = %d -- err= %10.6e\n",N,s); 
//-------------------- Free Host arrays
 free(A); 
 free(B);
 free(C);
//-------------------- Free GPU memory
 cudaFree(Ad);
 cudaFree(Bd);
 cudaFree(Cd);	
}

//-------------------- Prints error error Msg and exits 
void err_exit(char *errMsg) {
	printf("%s\n", errMsg);
	exit(1);
}
