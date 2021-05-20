/*
 * This source code is distributed under the terms defined  
 * in the file bfstdcooc_main.c of this source distribution.
*/
/* 
*  Breadth first search (BFS) 
*  Single precision (float data type) 
*  TurboBFS_COOC_TD:bfstdcooc.h
* 
*   This program computes: 
*   1) a sequential top-down BFS and
*   2) a single precision GPU-based parallel top-down BFS for 
*      unweighted graphs represented by sparse matrices in 
*      the COOC format.
*      
* 
*/

#ifndef BFSTDCOOC_H
#define BFSTDCOOC_H

#define MAX_THREADS_PER_GRID (2**31)
#define THREADS_PER_WARP 32
#define THREADS_PER_BLOCK 1024
#define WARPS_PER_BLOCK (THREADS_PER_BLOCK/THREADS_PER_WARP)


/*define Structure of Arrays (SoA) for the sparse matrix A in the COOC format*/
struct Cooc{
      int   *ICOOC;
      int   *JCOOC;
};

/**************************************************************************/
/* 
 * function to compute a sequential top-down BFS for unweighted graphs,
 * represented by sparse adjacency matrices in COOC format.   
 *  
*/
int bfs_seq_td_cooc (int *I,int *J,float *sigma,int r,int nz,int n);

/**************************************************************************/
/* 
 * function to compute a gpu-based parallel top-down BFS (scalar) for unweighted graphs 
 * represented by sparse adjacency matrices in COOC format. 
 *  
 */
int  bfs_gpu_td_cooc_sc (int *I_h,int *J_h,float *sigma_h,int r,int nz,
			 int n,int repetition);

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
 * function to compute the multiplication of one vector times the complement 
 * of the other vector.
*/
int mult_vs (float *sigma,float *f,int n);

/**************************************************************************/
/* 
 * function to assign one vector (different than zero) to another vector. 
 * 
*/
int assign_v(float *f,float *f_t,int n); 

#endif
