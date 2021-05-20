/*
 * This source code is distributed under the terms defined  
 * in the file bfstdcooc_main.c of this source distribution.
*/
/* 
*  Breadth first search (BFS) 
*  Single precision (float data type) 
*  TurboBFS_COOC_TD:bfsgpug.cuh
* 
*  This program defines the prototypes of some functions 
*  used to compute  the GPU-based parallel  top-down BFS
*  for unweighted graphs represented by sparse 
*  adjacency matrices in the COOC format.
*
*/

/*************************prototype kernel*************************************/
__global__ void bfsFunctionsKernel (int *f_d,int *ft_d,float *sigma_d,
				    int *c,int n);

