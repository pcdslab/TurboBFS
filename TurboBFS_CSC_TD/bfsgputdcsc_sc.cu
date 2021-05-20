/*
 * This source code is distributed under the terms defined  
 * in the file bfstdcsc_main.c of this source distribution.
 */
/* 
 *  Breadth first search (BFS) 
 *  Single precision (float data type) 
 *  TurboBFS_CSC_TD:bfsgputdcsc_sc.cu
 * 
 *  This program computes the GPU-based parallel 
 *  top-down BFS (scalar) for unweighted graphs represented 
 *  by sparse adjacency matrices in the CSC format, including
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


/*************************prototype kernel*************************************/
__global__ void spMvUgCscScKernel (int *CP_d,int *IC_d,int *ft_d,int *f_d,
				   float *sigma_d,int j,int r,int n);
/******************************************************************************/

/* 
 * Function to compute a gpu-based parallel top-down BFS (scalar) for 
 * unweighted graphs represented by sparse adjacency matrices in CSC format,
 * including the computation of the S vector to store the depth at which each  
 * vertex is  discovered.
 *  
 */
int  bfs_gpu_td_csc_sc (int *IC_h,int *CP_h,int *S_h,float *sigma_h,int r,
			int nz,int n,int repetition){
  float t_spmv;
  float t_spmv_t = 0.0;
  float t_bfsfunctions;
  float t_bfsfunctions_t = 0.0;
  float t_sum = 0.0;
  float t_bfs_avg;
  int i,d = 0,dimGrid;
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
  for (i = 0; i<repetition; i++){
    *c = 1;
    d = 0;
    checkCudaErrors(cudaMemset(f_d,0,sizeof(*f_d)*n));
    checkCudaErrors(cudaMemset(sigma_d,0,sizeof(*sigma_d)*n));
    while (*c){
      d = d + 1;
      *c = 0;

      cudaEventRecord(start);
      spMvUgCscScKernel <<<dimGrid,THREADS_PER_BLOCK>>> (CP_d,IC_d,ft_d,f_d,sigma_d,d,r,n);
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
  printf("\nbfsgputdcsc_sc::d = %d,r = %d,t_sum=%lfms \n",d,r,t_sum);
  t_bfs_avg = t_sum/repetition;

  /*Copy device memory (sigma_d) to host memory (sigma_h)*/
  checkCudaErrors(cudaMemcpy(sigma_h,sigma_d, n*sizeof(*sigma_d),cudaMemcpyDeviceToHost));

  /*Copy device memory (S_d) to host memory (S_h)*/
  checkCudaErrors(cudaMemcpy(S_h,S_d, n*sizeof(*S_d),cudaMemcpyDeviceToHost));

  int print_t = 1;
  if (print_t){
    printf("bfsgputdcsc_sc::time f <-- fA d = %lfms \n",t_spmv_t/repetition);
    printf("bfsgputdcsc_sc::time time bfs functions d = %lfms \n", t_bfsfunctions_t/repetition);
    printf("bfsgputdcsc_sc::average time BFS d = %lfms \n",t_bfs_avg);
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
}//end bfs_gpu_td_csc_sc


/**************************************************************************/
/* 
 * if d = 1, initialize f(r) and sigma(r),
 * compute the gpu-based parallel sparse matrix-vector multiplication    
 * for sparse matrices in the CSC format, representing unweighted 
 * graphs. 
 */
__global__
void spMvUgCscScKernel (int *CP_d,int *IC_d,int *ft_d,int *f_d,
			float *sigma_d,int d,int r,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < n){
    //initialize f(r) and sigma(r)
    if (d == 1){
      f_d[r] = 1;
      sigma_d[r] = 1.0;
    }
    //compute spmv
    ft_d[i] = 0;
    if (sigma_d[i] < 0.01){
      int k;
      int start = CP_d[i];
      int end = CP_d[i+1];
      int sum = 0;
      for (k = start;k < end; k++){
	sum += f_d[IC_d[k]];
      }
      if (sum > 0.9) {
	ft_d[i] = sum;
      }
    }
  }
}//end spMvUgCscScKernel

/**************************************************************************/
/*
 * assign vector ft_d to vector f_d,
 * check that the vector f_d  has at least one nonzero element
 * add the vector f to vector sigma.
 * compute the S vector. 
 */
__global__
void bfsFunctionsKernel (int *f_d,int *ft_d,float *sigma_d,int *S_d,int *c,
			 int n,int d){

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < n){
    f_d[i] = 0;
    if (ft_d[i] > 0.9){
      *c = 1;
      f_d[i] = ft_d[i];
      sigma_d[i] += ft_d[i];
      S_d[i] = d;
    }
  }
}//end  bfsFunctionsKernel
