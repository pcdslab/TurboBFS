/*
 * This source code is distributed under the terms defined on 
 * in the file bfstdcooc_main.c of this source distribution.
*/
/* 
*   Breadth first search (BFS)  
*   Single precision (float data type)
*   TurboBFS_COOC_TD:utils_bfs.h
* 
*  Functions to:
*  1) print results
*  2) check the GPU computed values of BFS shortest paths by comparing
*     with the results of the sequential BFS.
* 
*/

#ifndef  UTILS_BFS_H
#define  UTILS_BFS_H


/* 
 * function to print results
*/
int  printBFS(int *IC,int *JC,float *sigma,float *sigma_hgpu,int nz, int n);

/******************************************************************************/		
/* 
 * function to check the GPU computed values of BFS shortest paths by comparing
 * with the results of the sequential BFS. 
*/
int bfs_check(float *sigma,float *sigma_hgpu,int n);

#endif
