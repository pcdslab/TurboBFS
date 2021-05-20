/* 
 *  Breadth first search (BFS) 
 *  Single precision (float data type) 
 * 
 *  This program computes some functions needed
 *  for the computation of the sequential top-down
 *  BFS
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spmvcsc_seq.h"
#include "bfsbucsc.h"

/* 
 * function to compute the sum of two vectors
 * 
 */
int sum_vs (float *sigma,float *f,int n){

  int k;

  for (k=0; k<n; k++){
    if (f[k] > 0.1){
      sigma[k] += f[k];
    }
  }

  return 0;
}//end sum_vs

/**************************************************************************/
/* 
 * function to check that a vector is 0
 * 
 */
int check_f (int *c,float *f,int n){

  int k;
  for (k=0; k<n; k++){
    if (f[k] > 0.9){
      *c = 1;
      return 0;
    }
  }

  return 0;
}//end check_f

/**************************************************************************/
/* 
 * function to compute the S vector to store the depth at which each vertex 
 * is discovered.
 * 
 */
int S_v (int *S,float *f,int n,int d){

  int k;
  for (k=0; k<n; k++){
    if (f[k] > 0.9){
      S[k] = d;
    }
  }

  return 0;
}//end S_v

/**************************************************************************/
/* 
 * function to compute the multiplication of one vector times the complement of the
 * other vector.
 */
int mult_vs (float *sigma,float *f,int n){

  int k;

  for (k=0; k<n; k++){
    if (sigma[k] > 0.1){
      f[k] = 0.0;
    }
  }

  return 0;
}//end  mult_vs

/**************************************************************************/
/* 
 * function to assign one vector (different than zero) to another vector 
 * 
 */
int assign_v(float *f,float *f_t,int n){

  int k;

  for (k=0; k<n; k++){
    if (f_t[k] > 0.1){
      f[k] = f_t[k];
    }else{
      f[k] =  0.0;
    }
  }

  return 0;
}//end assign_v
