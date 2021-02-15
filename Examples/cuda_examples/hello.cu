#include <stdio.h>
  __global__ void helloFromGPU(){
    printf("Hello World-Thread: %d\n",threadIdx.x);                                                      
  }

  int main(void) {
    helloFromGPU<<<1,16>>>();
    cudaDeviceSynchronize();
    return(0);
  }
