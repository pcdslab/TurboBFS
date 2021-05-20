/* 
 * 
 * This program computes the sequential sparse matrix-vector 
 * multiplication for undirected, unweighted graphs represented 
 * by sparse adjacency matrices in the CSC format.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spmvcsc_seq.h"

/* 
 * function to compute the sequential sparse matrix-vector multiplication for 
 * undirected, unweighted graphs represented by sparse adjacency matrices in
 * the CSC format.
 *   
 */
int spmv_seq_ug_csc (float *f,int *I,int *CP,float *f_t,int n){

  int i, k, start, end;
  float sum;

  for (i=0; i<n; i++){
    f_t[i] = 0.0;
    sum = 0.0;
    start = CP[i];
    end = CP[i+1];
    for (k=start; k<end; k++){
      sum += f[I[k]];
    }
    if (sum > 0.1){
      f_t[i] = sum;
    }
  }

  return 0;
}//end spmv_seq_ug_csc
