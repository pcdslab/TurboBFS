/*
 * This source code is distributed under the terms defined  
 * in the file bfstdcooc_main.c of this source distribution.
 */
/* 
 *  Breadth first search (BFS) 
 *  Single precision (float data type) 
 *  TurboBFS_COOC_TD:spmvcooc_seq.h
 * 
 * This program computes sequential sparse matrix-vector multiplication  
 * (y= A'x) for unweighted graphs, with the sparse matrices in the COOC 
 * format.
 *  
 *  
 */

#ifndef SPMVCOOC_SEQ_H
#define SPMVCOOC_SEQ_H

/* 
 * function to compute the sequential sparse matrix-vector multiplication  
 *  (y= A'x) for unweighted graphs, with the sparse matrices in the COOC 
 *  format.
*/
int spmv_seq_td_cooc (float *f,int *I,int *J,float *f_t,int nz);


#endif
