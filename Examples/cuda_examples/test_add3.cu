#define NITER 1
#include <stdio.h>
#include <cuda.h>
#include <sys/time.h>
	   
/*-------------------- CUDA kernel. does an add  */
   __global__ void vec_add3(float *a, float *b, int n) {
  int block_id = gridDim.x * blockIdx.y + blockIdx.x ; 
  int i_t = blockDim.x * block_id + threadIdx.x;
  if (i_t<n) 	       
      a[i_t] += b[i_t] ; 
 }    
/*-------------------- main will execute on host */

int main(void){
/* -------------------- host & device arrays */
  float *y_h, *y_d,  *x_h, *x_d;  
  cudaError_t err;
  float ferr;
/* -------------------- size of arrays */
  const int N = 140000;
  int block_size ;
  int num_blocks_x, num_blocks_y;
  int i,len = 128*1024;
  size_t size = N * sizeof(float);
/* -------------------- Allocate array on host */
  y_h = (float *)malloc(size);        
  x_h = (float *)malloc(size);        
/* -------------------- Allocate array on device */
  cudaMalloc((void **) &y_d, size);   
  cudaMalloc((void **) &x_d, size);   
/* -------------------- Initialize host array & 
                        copy it to device */  
  for (int i=0; i<len; i++) {
      y_h[i] = (float) i / (float)len;
      x_h[i] = (float) (len-i)/(float)len;
  }     
  cudaMemcpy(y_d,y_h,size,cudaMemcpyHostToDevice);
  cudaMemcpy(x_d,x_h,size,cudaMemcpyHostToDevice);
/* -------------------- invoke kernel        */
 block_size = 1024;
 num_blocks_x = 16;	
 num_blocks_y = 8;
 dim3 grid_size(num_blocks_x, num_blocks_y, 1);
 vec_add3 <<< grid_size,block_size >>> (x_d, y_d,len);
//-------------------- get last error message
    err = cudaGetLastError() ;
    if (err != cudaSuccess) 
    printf(" error:  %s\n",cudaGetErrorString(err));
/*-------------------- get result from device */
   cudaMemcpy(x_h,x_d,size,cudaMemcpyDeviceToHost);
/* --------------------/ Print 20 first results / */
   ferr = 0.0;
   for (i=0; i<len; i++) 
      ferr += (1.0 - y_h[i])*(1.0 - x_h[i]);
   printf(" Done len = %d  \n Abs. Error is : %f\n", len, ferr);

/* -------------------- free memory */ 
   free(x_h); cudaFree(x_d);
   free(y_h); cudaFree(y_d);
 }

