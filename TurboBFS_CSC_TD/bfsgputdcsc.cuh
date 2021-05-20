/*
 * This source code is distributed under the terms defined  
 * in the file bfstdcsc_main.c of this source distribution.
*/
/* 
*  Breadth first search (BFS) 
*  Single precision (float data type) 
*  TurboBFS_CSC_TD:bfsgputdcsc.cuh
* 
*  This program defines the prototypes of some 
*  functions  used to compute the GPU-based parallel 
*  top-down BFS for unweighted graphs represented  
*  by sparse adjacency matrices in the CSC format.
*
*/

/*************************prototype kernel*************************************/
__global__ void bfsFunctionsKernel (int *f_d,int *ft_d,float *sigma_d,int *S,
				    int *c,int n,int d);

