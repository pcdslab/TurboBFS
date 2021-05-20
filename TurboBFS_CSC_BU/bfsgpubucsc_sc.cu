/*
 * This source code is distributed under the terms defined  
 * in the file bfsbucsc_main.c of this source distribution.
 */
/* 
 *  Breadth first search (BFS) 
 *  Single precision (float data type) 
 *  TurboBFS_CSC_BU:bfsgpubucsc_sc.cu
 * 
 *  This program computes gpu-based parallel bottom-up 
 *  BFS (scalar) for  unweighted graphs represented by  
 *  sparse adjacency matrices in CSC format. This BFS
 *  computes the S vector to store the  depth at which  
 *  each vertex is discovered.  
 *
 */

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cmath>

//includes CUDA project
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

extern "C"{
                 #include "bfsbucsc.h"

}


/*************************prototype kernel*************************************/
__global__ void bfsBuScKernel (int *CP_d,int *I_d,int *S_d,int *c,int d,int n,
			       int r,int bbreak);
/******************************************************************************/

/* 
 * function to compute a gpu-based parallel bottom-up BFS (scalar) for 
 * unweighted graphs represented by sparse adjacency matrices in CSC format. 
 * This BFS computes the S vector to store the depth at which each vertex  
 * is discovered.
 *  
 */
int  bfs_gpu_bu_csc_sc (int *I_h,int *CP_h,int *S_h,int r,int nz,int n,
			int repetition, int bbreak){
  float t_bfs = 0.0;
  float t_bfs_t = 0.0;
  float t_bfs_avg = 0.0;
  int i,d,dimGrid;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /*Allocate device memory for the vector CP_d */
  cudaEventRecord(start);
  int *CP_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&CP_d),sizeof(*CP_d)*(n+1)));

  /*Copy host memory (CP_h) to device memory (CP_d)*/
  checkCudaErrors(cudaMemcpy(CP_d,CP_h,(n+1)*sizeof(*CP_d),cudaMemcpyHostToDevice));

  /*Allocate device memory for the vector I_d */
  int *I_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&I_d),sizeof(*I_d)*nz));
  /*Copy host memory (I_h) to device memory (I_d)*/
  checkCudaErrors(cudaMemcpy(I_d,I_h,nz*sizeof(*I_d),cudaMemcpyHostToDevice));

  /*Allocate device memory for the vector S_d, and set S_d to -1. */
  int *S_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&S_d),sizeof(*S_d)*n));

  /*allocate unified memory for integer variable c for control of while loop*/
  int *c;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&c),sizeof(*c)));

  /*computing BFS */
  dimGrid = (n + THREADS_PER_BLOCK)/THREADS_PER_BLOCK;
  for (i = 0; i<repetition; i++){
    *c = 1;
    d = 0;
    checkCudaErrors(cudaMemset(S_d,-1,sizeof(*S_d)*n));
    while (*c){     
      *c = 0;
      cudaEventRecord(start);
      bfsBuScKernel <<<dimGrid,THREADS_PER_BLOCK>>> (CP_d,I_d,S_d,c,d,n,r,bbreak);			       
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&t_bfs,start,stop);
      t_bfs_t += t_bfs;
      d++;
    }
  }
  printf("\nbfsgpubucsc_sc::d = %d,r = %d,t_bfs_t=%lfms \n",d,r,t_bfs_t);
  t_bfs_avg = t_bfs_t/repetition;

  /*Copy device memory (S_d) to host memory (S_h)*/
  checkCudaErrors(cudaMemcpy(S_h,S_d, n*sizeof(*S_d),cudaMemcpyDeviceToHost));

  int print_t = 1;
  if (print_t) printf("bfsgpubucsc_sc:average time BFS  = %lfms \n",t_bfs_avg);

  /*cleanup memory*/
  checkCudaErrors(cudaFree(CP_d));
  checkCudaErrors(cudaFree(I_d));
  checkCudaErrors(cudaFree(S_d));
  checkCudaErrors(cudaFree(c));
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));

  return 0;
}//end bfs_gpu_bu_csc_sc

/******************************************************************************/
/* 
 * function to compute the bottom-up BFS with a gpu-based parallel scalar 
 * algorithm  for sparse matrices in the CSC format, representing unweighted 
 * graphs. 
 */
__global__
void bfsBuScKernel (int *CP_d,int *I_d, int *S_d,int *c,int d,int n,int r,
                    int bbreak){
		    

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < n){
    //initialize S_d(r)
    if (d == 0) S_d[r] = d;
     
    //compute BFS
    if (S_d[i] == -1){
      int start = CP_d[i];
      int end = CP_d[i+1];
      int cont = 1;
      int k = start;
      while (cont && k < end){
	if (S_d[I_d[k]] == d){
	  S_d[i] = d+1;
	  *c = 1;
	  if (bbreak) cont = 0;
        }
	k++;
      }
    }
  }
}//end bfsBuScKernel

