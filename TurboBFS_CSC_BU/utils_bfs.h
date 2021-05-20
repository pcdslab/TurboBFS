/*
 * This source code is distributed under the terms defined on 
 * in the file bfsbucsc_main.c of this source distribution.
*/
/* 
*   Breadth first search (BFS)  
*   Single precision (float data type)
*   TurboBFS_CSC_BU:utils_bfs.h
* 
*  This program:
*  1) prints results
*  2) check the GPU computed values of the S vector by comparing
*     with the results of the sequential computed values of S.
* 
*/

#ifndef  UTILS_BFS_H
#define  UTILS_BFS_H

/* 
 * function to print results
*/
int  printBFS(int *I,int *CP,int *S,int *S_hgpu,float *sigma,int nz, int n);

/******************************************************************************/
/* 
 * function to check the GPU computed values of the S vector by comparing
 * with the results of the sequential computed values of S. 
*/
int S_check(int *S,int *S_hgpu,int n);
#endif
