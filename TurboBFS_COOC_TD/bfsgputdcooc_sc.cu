/*
 * This source code is distributed under the terms defined  
 * in the file bfstdcooc_main.c of this source distribution.
 */
/* 
 *  Breadth first search (BFS) 
 *  Single precision (float data type) 
 *  TurboBFS_COOC_TD:bfsgputdcooc_sc.cu
 * 
 *  This program computes a GPU-based parallel top-down BFS
 *  (scalar) for unweighted connected graphs represented 
 *  by sparse adjacency matrices in the COOC format.
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
#include "bfsgputdcooc.cuh"
extern "C"{
                 #include "bfstdcooc.h"

}


/*************************prototype kernel*************************************/
__global__ void spMvUgCoocScKernel (int *I_d,int *J_d,int *ft_d,int *f_d,
				   float *sigma_d,int d,int r,int nz);
/******************************************************************************/

/* 
 * function to compute a sequential top-down BFS for unweighted graphs,
 * represented by sparse adjacency matrices in COOC format.   
 *  
*/

int  bfs_gpu_td_cooc_sc (int *I_h,int *J_h,float *sigma_h,int r,int nz,int n,
			int repetition){
			
  float t_spmv;
  float t_spmv_t = 0.0;
  float t_bfsfunctions;
  float t_bfsfunctions_t = 0.0;
  float t_sum = 0.0;
  float t_bfs_avg;
  int i,d = 0,dimGrid,dimGridspmv;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

   /*Allocate device memory for the vector I_d */
  int *I_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&I_d),sizeof(*I_d)*nz));
  /*Copy host memory (I_h) to device memory (I_d)*/
  checkCudaErrors(cudaMemcpy(I_d,I_h,nz*sizeof(*I_d),cudaMemcpyHostToDevice));

  /*Allocate device memory for the vector J_d */
  int *J_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&J_d),sizeof(*J_d)*nz));
  /*Copy host memory (J_h) to device memory (J_d)*/
  checkCudaErrors(cudaMemcpy(J_d,J_h,nz*sizeof(*J_d),cudaMemcpyHostToDevice));

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
  dimGridspmv = (nz + THREADS_PER_BLOCK)/THREADS_PER_BLOCK;
  for (i = 0; i<repetition; i++){
    *c = 1;
    d = 0;
    checkCudaErrors(cudaMemset(f_d,0,sizeof(*f_d)*n));
    checkCudaErrors(cudaMemset(sigma_d,0,sizeof(*sigma_d)*n));
    while (*c){
      d = d + 1;
      *c = 0;

      cudaEventRecord(start);
      checkCudaErrors(cudaMemset(ft_d,0,sizeof(*ft_d)*n));
      spMvUgCoocScKernel <<<dimGridspmv,THREADS_PER_BLOCK>>> (I_d,J_d,ft_d,f_d,sigma_d,d,r,nz);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&t_spmv,start,stop);
      t_spmv_t += t_spmv;

      cudaEventRecord(start);
      bfsFunctionsKernel <<<dimGrid,THREADS_PER_BLOCK>>> (f_d,ft_d,sigma_d,c,n);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&t_bfsfunctions,start,stop);
      t_bfsfunctions_t += t_bfsfunctions;
      
      t_sum += t_spmv + t_bfsfunctions;
    }
  }
  printf("\nbfsgpugcooc_sc::bfs_gpu_ug_cooc_sc::d = %d,r = %d,t_sum=%lfms \n",d,r,t_sum);
  t_bfs_avg = t_sum/repetition;

  /*Copy device memory (sigma_d) to host memory (sigma_h)*/
  checkCudaErrors(cudaMemcpy(sigma_h,sigma_d, n*sizeof(*sigma_d),cudaMemcpyDeviceToHost));
 
  int print_t = 1;
  if (print_t){
    printf("bfsgpugtdcooc_sc::time f <-- A'f d = %lfms \n",t_spmv_t/repetition);
    printf("bfsgputdcooc_sc::time time bfs functions d = %lfms \n", t_bfsfunctions_t/repetition);
    printf("bfsgputdcooc_sc::average time BFS d = %lfms \n",t_bfs_avg);
  }

  /*cleanup memory*/
  checkCudaErrors(cudaFree(J_d));
  checkCudaErrors(cudaFree(I_d));
  checkCudaErrors(cudaFree(sigma_d));
  checkCudaErrors(cudaFree(f_d));
  checkCudaErrors(cudaFree(ft_d));
  checkCudaErrors(cudaFree(c));
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));

  return 0;
}//end bfs_gpu_td_csc_sc


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/* 
 * if d = 1, initialize f(r) and sigma(r),
 * compute the gpu-based parallel sparse matrix-vector multiplication    
 * for sparse matrices in the CSC format, representing unweighted 
 * graphs. 
 */
 
__global__
void spMvUgCoocScKernel (int *I_d,int *J_d,int *ft_d,int *f_d,
			  float *sigma_d,int d,int r,int nz){

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < nz){
    if (d == 1){
      f_d[r] = 1;
      sigma_d[r] = 1.0;
    }
    if (f_d[I_d[i]] != 0) {
      int f = f_d[I_d[i]];
      atomicAdd(&ft_d[J_d[i]],f);
    }
  }
}//end spMvUgCoocScKernel

/******************************************************************************/
/*
 * if sigma(i) == 0, assign vector ft_d to vector f_d,
 * check that the vector f_d  has at least one nonzero element
 * add the vector f to vector sigma.
 */
__global__
void bfsFunctionsKernel (int *f_d,int *ft_d,float *sigma_d,int *c,
			 int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < n){
    f_d[i] = 0;
    if (sigma_d[i] < 0.01) f_d[i] = ft_d[i];
    if (f_d[i] > 0){
      *c = 1;      
      sigma_d[i] += f_d[i];
    }
  }
}//end  bfsFunctionsKernel
