/* 
*  Breadth first search (BFS) 
*  Single precision (float data type) 
* 
*  This program computes a sequential top-down
*  BFS for unweighted graphs represented 
*  by sparse matrices in the CSC format.
*
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "timer.h"
#include "spmvcsc_seq.h"
#include "bfstdcsc.h"

/* 
 * Function to compute a sequential top-down BFS for unweighted graphs,
 * represented by sparse adjacency matrices in CSC format, including the  
 * computation of the S vector to store the depth at which each vertex is
 * discovered.    
 *  
*/

int bfs_seq_td_csc (int *IC,int *CP,int *S,float *sigma,int r,int nz,int n){

    int d = 0;
    int c = 1;
    float *f;
    float *f_t;
    f =  (float *) calloc(n,sizeof(*f));
    f_t =  (float *) calloc(n,sizeof(*f_t));
    f[r] = 1;

    /* timing variables  */
    double initial_t;
    double delta;
    double sum_vs_t = 0;
    double spmv_t = 0;
    double assign_v_t = 0;
    double mult_vs_t = 0;
    double S_v_t = 0;
    double check_f_t = 0;
    double total_t;
    
    while (c > 0) {
      d++;
      initial_t = get_time();
      sum_vs (sigma,f,n);
      delta = get_time()-initial_t;
      sum_vs_t += delta;
      
      initial_t = get_time();
      spmv_seq_ug_csc (f,IC,CP,f_t,n);
      delta = get_time()-initial_t;
      spmv_t += delta;
      
      initial_t = get_time();
      assign_v(f,f_t,n);
      delta = get_time()-initial_t;
      assign_v_t += delta;
      
      initial_t = get_time();
      mult_vs (sigma,f,n);
      delta = get_time()-initial_t;
      mult_vs_t += delta;

      initial_t = get_time();
      S_v (S,f,n,d);
      delta = get_time()-initial_t;
      S_v_t += delta;
      
      initial_t = get_time();
      c = 0;
      check_f(&c,f,n);
      delta = get_time()-initial_t;
      check_f_t += delta;
    }

    printf("\nbfs_seq_ug_csc::d = %d,r = %d \n",d,r);
    total_t =  sum_vs_t +  spmv_t + assign_v_t +  mult_vs_t + S_v_t + check_f_t; 
    int p_t = 1;
    if (p_t) {    
      printf("bfstdcsc_seq::f <-- f +sigma time = %lfs \n",sum_vs_t);
      printf("bfstdcsc_seq::f_t <-- fA time = %lfs \n",spmv_t);
      printf("bfstdcsc_seq::f <-- f_t time = %lfs \n",assign_v_t);
      printf("bfstdcsc_seq::f <-- f*(-sigma) time = %lfs \n",mult_vs_t);
      printf("bfstdcsc_seq::S vector time = %lfs \n",S_v_t);
      printf("bfstdcsc_seq::c <-- check (f=0)) time = %lfs \n",check_f_t);
      printf("bfstdcsc_seq::sequential top-down BFS total time = %lfs \n",total_t);
    }
    
    return 0;
}//end  bfs_seq_td_csc
