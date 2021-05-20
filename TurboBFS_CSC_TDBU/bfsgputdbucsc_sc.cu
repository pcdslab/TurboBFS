/*
 * This source code is distributed under the terms defined  
 * in the file bfstdbucsc_main.c of this source distribution.
 */
/* 
 *  Breadth first search (BFS) 
 *  Single precision (float data type) 
 *  TurboBFS_CSC_TDBU:bfsgputdbucsc_sc.cu
 * 
 *  This program computes the GPU-based parallel combined 
 *  top-down/bottom-up BFS (scalar) for unweighted graphs
 *  represented  by sparse adjacency matrices in the CSC
 *  format.This combined BFS is a direction-optimization 
 *  algorithm switching from the top-down BFS to the bottom-up
 *  BFS when the frontier vector becomes dense.  
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
                     #include "bfstdbucsc.h"

}


/*************************prototype kernel*************************************/
__global__ void BfsCscScKernel (int *CP_d,int *I_d,int *ft_d,int *f_d,
				int *sigma_d,int *S_d,int d,int *c,int r,int n);
/******************************************************************************/
__global__ void fKernel (int *f_d,int *ft_d,float *sizef,int n);
/******************************************************************************/
__global__ void bfsBuScKernel (int *CP_d,int *I_d,int *S_d,int *c,int d,int n,
			       int r,int bbreak);
/******************************************************************************/

/* 
 * function to compute a gpu-based parallel BFS (scalar) for unweighted graphs 
 * represented by sparse adjacency matrices in CSC format.The BFS is a
 * direction-optimization algorithm switching from the top-down BFS to the
 * bottom-up BFS when the frontier vector becomes dense.
 *  
 */
int  bfs_gpu_tdbu_csc_sc (int *I_h,int *CP_h,int *S_h,int *sigma_h,int r,
			  int nz,int n,int repetition,int bbreak){
  float t_bfsk = 0.0;
  float t_bfsk_t = 0.0;
  float t_fk = 0.0;
  float t_fk_t = 0.0;
  float t_bfsbu = 0.0;
  float t_bfsbu_t = 0.0;
  float t_sum = 0.0;
  float t_bfs_avg = 0.0;
  float thr_td,thresh_td;
  int i,d = 0,dimGrid,count_td  = 0,count_bu  = 0;
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

  /*Allocate device memory for the vector S_d, and set S_d to zero. */
  int *S_d;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&S_d),sizeof(*S_d)*n));
  checkCudaErrors(cudaMemset(S_d,-1,sizeof(*S_d)*n));

  /*Allocate device memory for the vector sigma_d */
  int *sigma_d;
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

  /*allocate unified memory for integer variable sizef for storing the size of f*/
  float *sizef;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&sizef),sizeof(*sizef)));

  /*computing BFS */
  dimGrid = (n + THREADS_PER_BLOCK)/THREADS_PER_BLOCK;
  for (i = 0; i<repetition; i++){
    *c = 1;
    d = 0;
    *sizef = 0.0;
    thr_td = 0.001; //set by the user
    thresh_td = thr_td*n;
    count_td = 0;
    count_bu = 0;
    checkCudaErrors(cudaMemset(f_d,0,sizeof(*f_d)*n));
    checkCudaErrors(cudaMemset(sigma_d,0,sizeof(*sigma_d)*n));
    checkCudaErrors(cudaMemset(S_d,-1,sizeof(*S_d)*n));
    while (*c){

      *c = 0;
      if( *sizef <= thresh_td){
        *sizef = 0.0;
	d++;
	t_bfsk = 0.0;
	cudaEventRecord(start);
	BfsCscScKernel <<<dimGrid,THREADS_PER_BLOCK>>> (CP_d,I_d,ft_d,f_d,sigma_d,S_d,d,c,r,n);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&t_bfsk,start,stop);
	t_bfsk_t += t_bfsk;

	cudaEventRecord(start);
	fKernel <<<dimGrid,THREADS_PER_BLOCK>>> (f_d,ft_d,sizef,n);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&t_fk,start,stop);
	t_fk_t += t_fk;
	t_sum += t_bfsk + t_fk;
	count_td++;
      }else{
	t_bfsbu = 0.0;
	cudaEventRecord(start);
	bfsBuScKernel <<<dimGrid,THREADS_PER_BLOCK>>> (CP_d,I_d,S_d,c,d,n,r,bbreak);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&t_bfsbu,start,stop);
	t_bfsbu_t += t_bfsbu;
	t_sum += t_bfsbu;
	d++;
	count_bu++;
      }
    }
  }
  printf("\nbfsgputdbucsc_sc::d=%d,r=%d,thr_td=%lf,count_td=%d,count_bu=%d,t_sum=%lfms \n",d,r,thr_td,count_td,count_bu,t_sum);
  t_bfs_avg = t_sum/repetition;

  /*Copy device memory (sigma_d) to host memory (sigma_h)*/
  checkCudaErrors(cudaMemcpy(sigma_h,sigma_d, n*sizeof(*sigma_d),cudaMemcpyDeviceToHost));

  /*Copy device memory (S_d) to host memory (S_h)*/
  checkCudaErrors(cudaMemcpy(S_h,S_d, n*sizeof(*S_d),cudaMemcpyDeviceToHost));

  int print_t = 1;
  if (print_t){
    printf("bfsgputdbucsc_sc::time time bfs kernel = %lfms \n", t_bfsk_t/repetition);
    printf("bfsgputdbucsc_sc::time time f <-- ft kernel = %lfms \n", t_fk_t/repetition);
    printf("bfsgputdbucsc_sc::time time BFS TD  = %lfms \n", (t_bfsk_t+t_fk_t)/repetition);
    printf("bfsgputdbucsc_sc::time time BFS BU  = %lfms \n", t_bfsbu_t/repetition);
    printf("bfsgputdbucsc_sc::average time BFS = %lfms \n",t_bfs_avg);
  }

  /*cleanup memory*/
  checkCudaErrors(cudaFree(CP_d));
  checkCudaErrors(cudaFree(I_d));
  checkCudaErrors(cudaFree(S_d));
  checkCudaErrors(cudaFree(sigma_d));
  checkCudaErrors(cudaFree(f_d));
  checkCudaErrors(cudaFree(ft_d));
  checkCudaErrors(cudaFree(c));
  checkCudaErrors(cudaFree(sizef));
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));

  return 0;
}//end bfs_gpu_tdbu_csc_sc


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/* 
 * Function to:
 * 1)if d = 1, initialize f(r) and sigma(r),
 * 2)compute the gpu-based parallel vector-sparse matrix multiplication    
 * for sparse matrices in the CSC format, representing unweighted 
 * graphs.
 * 3) Update S, sigma vectors and c variable. 
 */
__global__
void BfsCscScKernel (int *CP_d,int *IC_d,int *ft_d,int *f_d,int *sigma_d,
		     int *S_d,int d,int *c,int r,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < n){
    //initialize f(r) and sigma(r)
    if (d == 1){
      f_d[r] = 1;
      sigma_d[r] = 1;
      S_d[r] = 0;
    }
    //compute bfs vectors
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
	sigma_d[i] += ft_d[i];
	if(S_d[i] == -1) S_d[i] = d;
	*c = 1;
      }
    }
  }
}//end BfsCscScKernel

/******************************************************************************/
/*
 * Function to assign vector ft_d != 0, to vector f_d,
 */
__global__
void fKernel (int *f_d,int *ft_d,float *sizef,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < n){
    f_d[i] = 0;
    if (ft_d[i] > 0.9){
      f_d[i] = ft_d[i];
      atomicAdd(sizef,1.0f);
    }
  }
}//end fKernel

////////////////////////////////////////////////////////////////////////////////
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
    if (d == 0) S_d[r] = d;

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
