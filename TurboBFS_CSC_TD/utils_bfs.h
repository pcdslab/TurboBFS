/*
 * This source code is distributed under the terms defined on 
 * in the file bfstdcsc_main.c of this source distribution.
*/
/* 
*   Breadth first search (BFS)  
*   Single precision (float data type)
*   TurboBFS_CSC_TD:utils_bfs.h
* 
*  Functions to:
*  1) print results
*  2) check the GPU computed values of BFS shortest paths by comparing
*     with the results of the sequential BFS.
*  3) check the GPU computed values of the S vector by comparing
*     with the results of the sequential computed values of S.
* 
*/

#ifndef  UTILS_BFS_H
#define  UTILS_BFS_H

/* 
 * function to print results
*/
int  printBFS(int *IC,int *CP,int *S,int *S_hgpu,float *sigma,float *sigma_hgpu,
	      int nz, int n);
	      
/******************************************************************************/		
/* 
 * function to check the GPU computed values of BFS shortest paths by comparing
 * with the results of the sequential BFS. 
*/
int bfs_check(float *sigma,float *sigma_hgpu,int n);

/******************************************************************************/
/* 
 * function to check the GPU computed values of the S vector by comparing
 * with the results of the sequential computed values of S. 
*/
int S_check(int *S,int *S_hgpu,int n);
#endif
