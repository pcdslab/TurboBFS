
/* 
 * 
 * This program computes the sequential sparse matrix-vector 
 * multiplication for undirected, unweighted graphs represented 
 * by sparse adjacency matrices in the CSC format.
 * integer data type 
 * 
 */

#ifndef SPMVCSC_SEQ_H
#define SPMVCSC_SEQ_H

/* 
 * function to compute the sequential sparse matrix-vector  multiplication for 
 * undirected, unweighted graphs represented by sparse adjacency matrices in
 * the CSC format.
 *   
 */
int spmv_seq_ug_csc (int *f,int *I,int *CP,int *f_t,int n);


#endif
