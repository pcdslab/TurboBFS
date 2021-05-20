/*
 * This source code is distributed under the terms defined  
 * in the file bfstdbucsc_main.c of this source distribution.
 */
/* 
 *  Breadth first search (BFS) 
 *  Single precision (float data type) 
 *  TurboBFS_CSC_TDBU:bfsgputdbucsc_wa.cu
 * 
 *  This program computes the GPU-based parallel combined 
 *  top-down/bottom-up BFS ((warp shuffle) for unweighted 
 *  graphs represented by sparse adjacency matrices in the 
 *  CSC format.The BFS is a direction-optimization algorithm 
 *  switching from the  top-down BFS to the bottom-up BFS  
 *  when the frontier vector becomes dense.
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


/*************************prototype kernel MVS*********************************/
__global__ void spMvUgCscWaKernel (int *CP_d,int *I_d,int *ft_d,int *f_d,
				   int *sigma_d,int *S_d,int d,int r,int n);
/******************************************************************************/
__global__ void bfsBuWaKernel(int *CP_d,int *I_d, int *S_d,int *c,int d,
			      int n,int r,int bbreak);
/******************************************************************************/

__global__ void bfsTdBuFunctionsKernel (int *f_d,int *ft_d,int *sigma_d,int *S_d,
				        int *c,float *sizef,int n,int d);
/******************************************************************************/

/* 
 * function to compute a gpu-based parallel BFS (warp shuffle) for unweighted
 * graphs represented by sparse adjacency matrices in CSC format.The BFS is a
 * direction-optimization algorithm switching from the top-down BFS to the
 * bottom-up BFS when the frontier vector becomes dense.
 *  
 */
int bfs_gpu_tdbu_csc_wa (int *I_h,int *CP_h,int *S_h,int *sigma_h,int r,
			 int nz,int n,int repetition,int bbreak){
  float t_spmv = 0.0;
  float t_spmv_t = 0.0;
  float t_bfsfunctions = 0.0;
  float t_bfsfunctions_t = 0.0;
  float t_bfsbu = 0.0;
  float t_bfsbu_t = 0.0;
  float t_sum = 0.0;
  float t_bfs_avg = 0.0;
  float thr_td,thresh_td;
  int i,d,dimGrid,count_td  = 0,count_bu = 0;
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
  dimGrid = (n + THREADS_PER_WARP)/THREADS_PER_WARP;
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
	t_spmv = 0.0;
	cudaEventRecord(start);
	checkCudaErrors(cudaMemset(ft_d,0,sizeof(*ft_d)*n));
	spMvUgCscWaKernel <<<dimGrid,THREADS_PER_BLOCK>>> (CP_d,I_d,ft_d,f_d,sigma_d,S_d,d,r,n);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&t_spmv,start,stop);
	t_spmv_t += t_spmv;

	t_bfsfunctions = 0.0;
	cudaEventRecord(start);
	bfsTdBuFunctionsKernel <<<dimGrid,THREADS_PER_BLOCK>>> (f_d,ft_d,sigma_d,S_d,c,sizef,n,d);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&t_bfsfunctions,start,stop);
	t_bfsfunctions_t += t_bfsfunctions;
	t_sum += t_spmv + t_bfsfunctions;
	count_td++;
      }else{
      	t_bfsbu = 0.0;
	cudaEventRecord(start);
	bfsBuWaKernel <<<dimGrid,THREADS_PER_BLOCK>>> (CP_d,I_d,S_d,c,d,n,r,bbreak);
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
  printf("\nbfsgputdbucsc_wa::d=%d,r=%d,thr_td=%lf,count_td=%d,count_bu=%d,t_sum=%lfms \n",d,r,thr_td,count_td,count_bu,t_sum);
  t_bfs_avg = t_sum/repetition;

  /*Copy device memory (sigma_d) to host memory (sigma_h)*/
  checkCudaErrors(cudaMemcpy(sigma_h,sigma_d, n*sizeof(*sigma_d),cudaMemcpyDeviceToHost));

  /*Copy device memory (S_d) to host memory (S_h)*/
  checkCudaErrors(cudaMemcpy(S_h,S_d, n*sizeof(*S_d),cudaMemcpyDeviceToHost));

  int print_t = 1;
  if (print_t){
    printf("\nbfsgputdbucsc_wa::time f <-- fA = %lfms \n",t_spmv_t/repetition);
    printf("bfsgputdbucsc_wa::time time bfs functions = %lfms \n", t_bfsfunctions_t/repetition);
    printf("bfsgputdbucsc_wa::time time BFS TD  = %lfms \n", (t_bfsfunctions_t+t_spmv_t)/repetition);
    printf("bfsgputdbucsc_wa::time time BFS BU  = %lfms \n", t_bfsbu_t/repetition);
    printf("bfsgputdbucsc_wa::average time BFS  = %lfms \n",t_bfs_avg);
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
  //checkCudaErrors(cudaFree(sizeS));
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));

  return 0;
}//end bfs_gpu_tdbu_csc_wa


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/* 
 * if d = 1, initialize f(r) and sigma(r),
 * compute the gpu-based parallel sparse matrix-vector multiplication    
 * for sparse matrices in the CSC format, representing unweighted 
 * graphs. 
 */
__global__
void spMvUgCscWaKernel (int *CP_d,int *IC_d,int *ft_d,int *f_d,
			int *sigma_d,int *S_d,int d,int r,int n){

  __shared__ volatile int cp [WARPS_PER_BLOCK][2];


  //initialize f(r),sigma(r),S(r)
  if (d == 1){
    f_d[r] = 1;
    sigma_d[r] = 1;
    S_d[r] = 0;
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*
 * assign vector ft_d to vector f_d,
 * check that the vector f_d  has at least one nonzero element
 * add the vector f to vector sigma.
 * computes the S vector
 * computes the size of f 
 */
__global__
void bfsTdBuFunctionsKernel (int *f_d,int *ft_d,int *sigma_d,int *S_d,
			   int *c,float *sizef,int n,int d){

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < n){
    f_d[i] = 0;
    if (ft_d[i] > 0.9){
      *c = 1;
      atomicAdd(sizef,1.0f);
      f_d[i] = ft_d[i];
      sigma_d[i] += ft_d[i];
      if (S_d[i] == -1) S_d[i] = d;
    }
  }

}//end  bfsFunctionsKernel

////////////////////////////////////////////////////////////////////////////////
/* 
 * function to compute the bottom-up BFS with a gpu-based parallel (warp-based)
 * algorithm for sparse matrices in the CSC format, representing unweighted 
 * graphs. 
 */
__global__
void bfsBuWaKernel (int *CP_d,int *I_d, int *S_d,int *c,int d,
		    int n,int r,int bbreak){

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
  int cont = 1;
  
  while (cont && col < n  && S_d[col] == -1) {
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
	  if (bbreak) cont = 0;
	}
      }
      icp += THREADS_PER_WARP;
      while (cont && icp < end){
	if(S_d[I_d[icp]] == d){
	  S_d[col] = d+1;
	  *c = 1;
	  if (bbreak)  cont = 0;
	}
	icp += THREADS_PER_WARP;
      }
    }else{//number of column elements <= THREADS_PER_WARP){
      icp = start + thread_lane_id;
      while (cont && icp < end){
	if(S_d[I_d[icp]] == d){
	  S_d[col] = d+1;
	  *c = 1;
	  if (bbreak) cont = 0;
	}
	icp += THREADS_PER_WARP;
      }
    }
    col += num_warps;
  }
}//bfsBuWaKernel
