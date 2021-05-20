/*
 * This source code is distributed under the terms defined  
 * in the file bfsbucsc_main.c of this source distribution.
*/
/* 
 *  Breadth first search (BFS) 
 *  Single precision (float data type) 
 *   TurboBFS_CSC_BU:bfsgpubucsc_wa.cu
 * 
 *  This program computes the GPU-based parallel parallel 
 *  Bottom-Up BFS (warp shuffle) for  unweighted graphs
 *  represented by sparse adjacency matrices in CSC format. 
 *  This BFS computes the S vector to store the depth at 
 *  which each vertex is discovered.
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


/*************************prototype kernel MVS*********************************/
__global__ void bfsBuWaKernel(int *CP_d,int *I_d, int *S_d,int *c,int d,int n,
	   	               int r,int bbreak);
/******************************************************************************/

/* 
 * function to compute a gpu-based parallel bottom-up BFS (warp shuffle) for 
 * unweighted graphs represented by sparse adjacency matrices in CSC format. 
 * This BFS computes the S vector to store the depth at which each vertex  
 * is discovered.
 *  
 */
int  bfs_gpu_bu_csc_wa (int *I_h,int *CP_h,int *S_h,int r,int nz,int n,
			int repetition,int bbreak){
  float t_bfs = 0.0;
  float t_bfs_t = 0.0;
  float t_bfs_avg = 0.0;
  int i,d,dimGrid;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /*Allocate device memory for the vector CP_d */
  int *CP_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&CP_d),sizeof(*CP_d)*(n+1)));
  
  /*Copy host memory (CP_h) to device memory (CP_d)*/
  checkCudaErrors(cudaMemcpy(CP_d,CP_h,(n+1)*sizeof(*CP_d),cudaMemcpyHostToDevice));

  /*Allocate device memory for the vector I_d */
  int *I_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&I_d),sizeof(*I_d)*nz));
  /*Copy host memory (I_h) to device memory (I_d)*/
  checkCudaErrors(cudaMemcpy(I_d,I_h,nz*sizeof(*I_d),cudaMemcpyHostToDevice));

  /*Allocate device memory for the vector S_d */
  int *S_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&S_d),sizeof(*S_d)*n));
 
  /*allocate unified memory for integer variable c for control of while loop*/
  int *c;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&c),sizeof(*c)));
  
  /*computing BFS */
  dimGrid = (n + THREADS_PER_WARP)/THREADS_PER_WARP;
  for (i = 0; i<repetition; i++){
    *c = 1;
    d = 0;
    checkCudaErrors(cudaMemset(S_d,-1,sizeof(*S_d)*n));
    while (*c){      
      *c = 0;
      
      cudaEventRecord(start);
      bfsBuWaKernel <<<dimGrid,THREADS_PER_BLOCK>>> (CP_d,I_d,S_d,c,d,n,r,bbreak);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&t_bfs,start,stop);
      t_bfs_t += t_bfs;
      d = d+1;
    }
  }
  printf("\nbfsgpubucsc_wa::d = %d,r = %d,t_bfs_t=%lfms \n",d,r,t_bfs_t);
  t_bfs_avg = t_bfs_t/repetition;

  /*Copy device memory (S_d) to host memory (S_h)*/
  checkCudaErrors(cudaMemcpy(S_h,S_d, n*sizeof(*S_d),cudaMemcpyDeviceToHost));

  int print_t = 1;
  if (print_t) printf("bfsgpubucsc_wa:average time BFS  = %lfms \n",t_bfs_avg);

  /*cleanup memory*/
  checkCudaErrors(cudaFree(CP_d));
  checkCudaErrors(cudaFree(I_d));
  checkCudaErrors(cudaFree(S_d));
  checkCudaErrors(cudaFree(c));
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));

  return 0;
}//end bfs_gpu_bu_csc_wa

/******************************************************************************/
/* 
 * function to compute the Bottom-Up BFS with a gpu-based parallel (warp-based)
 * algorithm for sparse matrices in the CSC format, representing unweighted 
 * graphs. 
 */
__global__
void bfsBuWaKernel (int *CP_d,int *I_d,int *S_d,int *c,int d,int n,int r,
		    int bbreak){			       
				   
  __shared__ volatile int cp [WARPS_PER_BLOCK][2];
  
  //initialize S_d(r)
  if (d == 0) S_d[r] = d;

  //compute bfs
  int thread_id = threadIdx.x + blockIdx.x * blockDim.x; //global thread index
  int thread_lane_id = thread_id & (THREADS_PER_WARP-1); //thread index within the warp
  int warp_id = thread_id/THREADS_PER_WARP; //global warp index
  int warp_lane_id = threadIdx.x/ THREADS_PER_WARP; //warp index within the block
  int num_warps =  WARPS_PER_BLOCK*gridDim.x;       //total number of available warps

  int col = warp_id;
  int icp;
  int brk = 1;
  
  while (brk && col < n  && S_d[col] == -1) { 
    if (thread_lane_id<2){
      cp[warp_lane_id][thread_lane_id] = CP_d[col+thread_lane_id];
    }
    int start = cp[warp_lane_id][0];
    int end = cp[warp_lane_id][1];
      
    //compute S
    if (end - start > THREADS_PER_WARP){//number of column elements > THREADS_PER_WARP
      icp = start -(start & (THREADS_PER_WARP-1)) + thread_lane_id;
      if (icp >= start && icp < end){
	if(S_d[I_d[icp]] == d){
	  S_d[col] = d+1;
	  *c = 1;
	  if (bbreak) brk = 0;
	}
      }
      icp += THREADS_PER_WARP;
      while (brk && icp < end){
        if(S_d[I_d[icp]] == d){
	  S_d[col] = d+1;
	  *c = 1;
	  if (bbreak)  brk = 0;
	}
	icp += THREADS_PER_WARP;
      }
      }else{//number of column elements <= THREADS_PER_WARP
	icp = start + thread_lane_id;
	while (brk && icp < end){
	  if(S_d[I_d[icp]] == d){
	    S_d[col] = d+1;
	    *c = 1;
	    if (bbreak) brk = 0;
	  }
	  icp += THREADS_PER_WARP;
	}
      }
    col += num_warps;
  }
}//bfsBuWaKernel
