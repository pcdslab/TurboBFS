/*
 * This source code is distributed under the terms defined  
 * in the file bfstdcooc_main.c of this source distribution.
 */
/* 
 *  Breadth first search (BFS) 
 *  Single precision (float data type) 
 *  TurboBFS_COOC_TD:spmvcooc_seq.c
 * 
 * This program computes:sequential sparse matrix-vector multiplication  
 * (y= A'x) for unweighted graphs, with the sparse matrices in the COOC 
 * format.
 *  
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spmvcooc_seq.h"

/* 
 * function to compute the sequential sparse matrix-vector multiplication  
 *  (y= A'x) for unweighted graphs, with the sparse matrices in the COOC 
 *  format.
*/
int spmv_seq_td_cooc (float *f,int *I,int *J,float *f_t,int nz){

  int i;

    for (i=0; i<nz; i++){
       f_t[J[i]] += f[I[i]];
    }

  return 0;
}//end spmv_seq_td_csc
