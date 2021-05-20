/*
 * This source code is distributed under the terms defined  
 * in the file bfsbucsc_main.c of this source distribution.
*/
/* 
*  Breadth first search (BFS) 
*  Single precision (float data type) 
*  TurboBFS_CSC_BU:bfsbucsc.h
* 
*   This program computes: 
*   1) a sequential top-down BFS and
*   2) a single precision GPU-based parallel bottom-up BFS 
*      for unweighted graphs represented with sparse adjacency  
*      matrices in the CSC format.
*      
* 
*/

#ifndef BFSBUCSC_H
#define BFSBUCSC_H

#define MAX_THREADS_PER_GRID (2**31)
#define THREADS_PER_WARP 32
#define THREADS_PER_BLOCK 1024
#define WARPS_PER_BLOCK (THREADS_PER_BLOCK/THREADS_PER_WARP)
#define I_SIZE ((3/2)*THREADS_PER_BLOCK)


/*define Structure of Arrays (SoA) for the sparse matrix A representing
 unweighted graphs in the CSC format*/
struct Csc{
  int   *IC;
  int   *CP;
  
};

/**************************************************************************/
/* 
 * function to compute a sequential top-down BFS for unweighted graphs,
 * represented by sparse adjacency matrices in CSC format.This includes 
 * the  computation of  the S vector to store the depth at which each    
 * vertex is discovered. 
 * 
*/
int bfs_seq_td_csc (int *I,int *CP,int *S,float *sigma,int r,int nz,int n);

/**************************************************************************/
/* 
 * function to compute a gpu-based parallel bottom-up BFS (scalar) for 
 * unweighted graphs represented by sparse adjacency matrices in CSC format. 
 * This BFS computes the S vector to store the depth at which each vertex  
 * is discovered.
 *  
 */
int  bfs_gpu_bu_csc_sc (int *I_h,int *CP_h,int *S_h,int r,int nz,int n,
			int repetition, int bbreak);

/**************************************************************************/
/* 
 * function to compute a gpu-based parallel bottom-up BFS (warp shuffle) for 
 * unweighted graphs represented by sparse adjacency matrices in CSC format. 
 * This BFS computes the S vector to store the depth at which each vertex  
 * is discovered.
 *  
 */
int  bfs_gpu_bu_csc_wa (int *I_h,int *CP_h,int *S_h,int r,int nz,int n,
			int repetition,int bbreak);

/**************************************************************************/
/* 
 * function to compute the sum of two vectors
 * 
*/
int sum_vs (float *sigma,float *f,int n);

/**************************************************************************/
/* 
 * function to check that a vector is 0
 * 
 */
int check_f (int *c,float *f,int n);

/**************************************************************************/
/* 
 * function to compute the S vector to store the depth at which each vertex 
 * is discovered.
 * 
 */
int S_v (int *S,float *f,int n,int d);

/**************************************************************************/
/* 
 * function to compute the multiplication of one vector times the complement of the
 * other vector.
*/
int mult_vs (float *sigma,float *f,int n);

/**************************************************************************/
/* 
 * function to assign one vector (different than zero) to another vector. 
 * 
*/
int assign_v(float *f,float *f_t,int n); 

#endif
