/*
 * This source code is distributed under the terms defined  
 * in the file bfstdbucsc_main.c of this source distribution.
*/
/* 
*  Breadth first search (BFS) 
*  Single precision (float data type) 
*  TurboBFS_CSC_TDBU:bfstdbucsc.h
* 
*   This program computes: 
*   1) a sequential top-down BFS and
*   2) computes the GPU-based parallel combined 
 *     top-down/bottom-up BFS for unweighted graphs
 *     represented  by sparse adjacency matrices in the CSC
 *     format.
*      
* 
*/

#ifndef BFSTDBUCSC_H
#define BFSTDBUCSC_H

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
 * function to compute a sequential top_down BFS for unweighted graphs,
 * represented by sparse adjacency matrices in CSC format, including
 * the  computation of the S vector to store the depth at which each 
 * vertex is discovered.    
 *  
*/
int bfs_seq_td_csc (int *IC,int *CP,int *S,int *sigma,int r,int nz,int n);

/**************************************************************************/
/* 
 * function to compute a gpu-based parallel BFS (scalar) for unweighted graphs 
 * represented by sparse adjacency matrices in CSC format.The BFS is a
 * direction-optimization algorithm switching from the top-down BFS to the
 * bottom-up BFS when the frontier vector becomes dense.
 *  
 */
int bfs_gpu_tdbu_csc_sc (int *I_h,int *CP_h,int *S_h,int *sigma_h,int r,
			 int nz,int n,int repetition,int bbreak);

/**************************************************************************/
/* 
 * function to compute a gpu-based parallel BFS (warp shuffle) for unweighted
 * graphs represented by sparse adjacency matrices in CSC format.The BFS is a
 * direction-optimization algorithm switching from the top-down BFS to the
 * bottom-up BFS when the frontier vector becomes dense.
 *  
 */
int bfs_gpu_tdbu_csc_wa (int *IC_h,int *CP_h,int *S_h,int *sigma_h,int r,
                         int nz,int n,int repetition,int bbreak);

/**************************************************************************/
/* 
 * function to compute the sum of two vectors
 * 
*/
int sum_vs (int *sigma,int *f,int n);

/**************************************************************************/
/* 
 * function to check that a vector is 0
 * 
 */
int check_f (int *c,int *f,int n);

/**************************************************************************/
/* 
 * function to compute the S vector to store the depth at which each vertex 
 * is discovered.
 * 
 */
int S_v (int *S,int *f,int n,int d);

/**************************************************************************/
/* 
 * function to compute the multiplication of one vector times the complement of the
 * other vector.
*/
int mult_vs (int *sigma,int *f,int n);

/**************************************************************************/
/* 
 * function to assign one vector (different than zero) to another vector. 
 * 
*/
int assign_v(int *f,int *f_t,int n); 

#endif
