/*
 * This source code is distributed under the terms defined  
 * in the file bfstdcsc_main.c of this source distribution.
*/
/* 
 *  Breadth first search (BFS) 
 *  Single precision (float data type) 
 *  TurboBFS_CSC_TD:bfsgputdcsc_wa.cu
 * 
 *  This program computes the GPU-based parallel top-down
 *  BFS (warp shuffle) for unweighted graphs represented by 
 *  sparse adjacency matrices in the CSC format, including
 *  the computation of the S array to store the depth at 
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
#include "bfsgputdcsc.cuh"

extern "C"{
             #include "bfstdcsc.h"
}

/*************************prototype kernel MVS*********************************/
__global__ void spMvUgCscWaKernel (int *CP_d,int *IC_d,int *ft_d,int *f_d,
				   float *sigma_d,int d,int r,int n);
/******************************************************************************/

/* 
 * function to compute a gpu-based parallel top-down BFS (warp shuffle) for
 * unweighted graphs represented by sparse adjacency matrices in CSC format,
 * including the computation of the S vector to store the depth at which each 
 * vertex is discovered.
 *  
 */
int  bfs_gpu_td_csc_wa (int *IC_h,int *CP_h,int *S_h,float *sigma_h,int r,
			int nz,int n,int repetition){
  float t_spmv;
  float t_spmv_t = 0.0;
  float t_bfsfunctions;
  float t_bfsfunctions_t = 0.0;
  float t_sum = 0.0;
  float t_bfs_avg;
  int i, d,  dimGrid_mvsp, dimGrid;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /*Allocate device memory for the vector CP_d */
  int *CP_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&CP_d),sizeof(*CP_d)*(n+1)));  
  /*Copy host memory (CP_h) to device memory (CP_d)*/
  checkCudaErrors(cudaMemcpy(CP_d,CP_h,(n+1)*sizeof(*CP_d),cudaMemcpyHostToDevice));

  /*Allocate device memory for the vector IC_d */
  int *IC_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&IC_d),sizeof(*IC_d)*nz));
  /*Copy host memory (IC_h) to device memory (IC_d)*/
  checkCudaErrors(cudaMemcpy(IC_d,IC_h,nz*sizeof(*IC_d),cudaMemcpyHostToDevice));

  /*Allocate device memory for the vector S_d, and set S_d to zero. */
  int *S_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&S_d),sizeof(*S_d)*n));
  checkCudaErrors(cudaMemset(S_d,0,sizeof(*S_d)*n));

  /*Allocate device memory for the vector sigma_d */
  float *sigma_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&sigma_d),sizeof(*sigma_d)*n));

  /*Allocate device memory for the vector f_d*/
  int *f_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&f_d),sizeof(*f_d)*n));

  /*Allocate device memory for the vector ft_d*/
  int *ft_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&ft_d),sizeof(*ft_d)*n));

  /*allocate unified memory for integer variable c for control of while loop*/
  int *c;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&c),sizeof(*c)));
  
  /*computing BFS */
  dimGrid = (n + THREADS_PER_BLOCK)/THREADS_PER_BLOCK;
  dimGrid_mvsp = (n + THREADS_PER_WARP)/THREADS_PER_WARP;
  for (i = 0; i<repetition; i++){
    *c = 1;
    d = 0;
    checkCudaErrors(cudaMemset(f_d,0,sizeof(*f_d)*n));
    checkCudaErrors(cudaMemset(sigma_d,0,sizeof(*sigma_d)*n));  
    while (*c){
      d = d+1;
      *c = 0;
      
      cudaEventRecord(start);
      checkCudaErrors(cudaMemset(ft_d,0,sizeof(*ft_d)*n));
      spMvUgCscWaKernel <<<dimGrid_mvsp,THREADS_PER_BLOCK>>> (CP_d,IC_d,ft_d,f_d,sigma_d,d,r,n);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&t_spmv,start,stop);
      t_spmv_t += t_spmv;
      
      cudaEventRecord(start);
      bfsFunctionsKernel <<<dimGrid,THREADS_PER_BLOCK>>> (f_d,ft_d,sigma_d,S_d,c,n,d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&t_bfsfunctions,start,stop);
      t_bfsfunctions_t += t_bfsfunctions;
      
      t_sum += t_spmv + t_bfsfunctions;
    }
  }
  printf("\nbfsgputdcsc_wa::d = %d,r = %d,t_sum=%lfms \n",d,r,t_sum);
  t_bfs_avg = t_sum/repetition;

  /*Copy device memory (sigma_d) to host memory (sigma_h)*/
  checkCudaErrors(cudaMemcpy(sigma_h,sigma_d, n*sizeof(*sigma_d),cudaMemcpyDeviceToHost));

  /*Copy device memory (S_d) to host memory (S_h)*/
  checkCudaErrors(cudaMemcpy(S_h,S_d, n*sizeof(*S_d),cudaMemcpyDeviceToHost));

  int print_t = 1;
  if (print_t){
    printf("bfsgputdcsc_wa::time f <-- fA = %lfms \n",t_spmv_t/repetition);
    printf("bfsgputdcsc_wa::time time bfs functions = %lfms \n", t_bfsfunctions_t/repetition);
    printf("bfsgputdcsc_wa::average time BFS  = %lfms \n",t_bfs_avg);
  }

  /*cleanup memory*/
  checkCudaErrors(cudaFree(CP_d));
  checkCudaErrors(cudaFree(IC_d));
  checkCudaErrors(cudaFree(S_d));
  checkCudaErrors(cudaFree(sigma_d));
  checkCudaErrors(cudaFree(f_d));
  checkCudaErrors(cudaFree(ft_d));
  checkCudaErrors(cudaFree(c));
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));

  return 0;
}//end bfs_gpu_td_csc_wa

/**************************************************************************/
/* 
 * if d = 1, initialize f(r) and sigma(r),
 * compute the gpu-based parallelsparse matrix-vector multiplication    
 * for sparse matrices in the CSC format, representing unweighted 
 * graphs. 
 */
__global__
void spMvUgCscWaKernel (int *CP_d,int *IC_d,int *ft_d,int *f_d,
			float *sigma_d,int d,int r,int n){
				   
  __shared__ volatile int cp [WARPS_PER_BLOCK][2];

  
  //initialize f(r) and sigma(r)
  if (d == 1){
      f_d[r] = 1;
      sigma_d[r] = 1.0;
  }
  //compute spmv
  int thread_id = threadIdx.x + blockIdx.x * blockDim.x; //global thread index
  int thread_lane_id = thread_id & (THREADS_PER_WARP-1); //thread index within the warp
  int warp_id = thread_id/THREADS_PER_WARP; //global warp index
  int warp_lane_id = threadIdx.x/ THREADS_PER_WARP; //warp index within the block
  int num_warps =  WARPS_PER_BLOCK*gridDim.x;       //total number of available warps

  int col;
  int icp;
  unsigned mask = 0xffffffff;

  for (col = warp_id; col < n; col += num_warps){
    if(sigma_d[col] < 0.01) {
    
      if (thread_lane_id<2){
        cp[warp_lane_id][thread_lane_id] = CP_d[col+thread_lane_id];
      }
      int start = cp[warp_lane_id][0];
      int end = cp[warp_lane_id][1];
      //printf("threadIdx.x=%d,blockIdx.x=%d,thread_lane_id=%d,warp_id=%d,warp_lane_id=%d,num_warps=%d,col=%d,start=%d,end=%d\n",threadIdx.x,blockIdx.x,thread_lane_id,warp_id,warp_lane_id,num_warps,col,start,end);
      int sum = 0;

      if (end - start > THREADS_PER_WARP){//number of column elements > THREADS_PER_WARP
        icp = start -(start & (THREADS_PER_WARP-1)) + thread_lane_id;
        /*accumulate local sums*/
        if (icp >= start && icp < end){
	  sum += f_d[IC_d[icp]];
        }
        /*accumulate local sums*/
        for (icp += THREADS_PER_WARP; icp < end; icp += THREADS_PER_WARP){
          sum += f_d[IC_d[icp]];
        }
      }else{//number of column elements <= THREADS_PER_WARP
        /*accumulate local sums*/
        for (icp = start + thread_lane_id; icp < end; icp += THREADS_PER_WARP){
          sum += f_d[IC_d[icp]];
        }
      }
      /*reduce local sums by the warp shuffle instruction */
      for (int offset = THREADS_PER_WARP/2; offset > 0; offset /= 2){
        sum += __shfl_down_sync(mask,sum,offset);
      }
      /*first thread in the warp output the final result*/
      if (thread_lane_id == 0){
        ft_d[col] = sum;
      }   
    }
  }
}//end spMvUgCscWaKernel
