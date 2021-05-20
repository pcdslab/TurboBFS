/*
 * This source code is distributed under the terms defined  
 * in the file bfstdcsc_main.c of this source distribution.
*/
/* 
*  Breadth first search (BFS) 
*  Single precision (float data type) 
*  TurboBFS_COOC_TD:bfstdcooc_seq.c
* 
*  This program computes a single precision top-down sequential 
*  BFS for unweighted graphs represented by sparse matrices
*  in the COOC format.
*
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "timer.h"
#include "spmvcooc_seq.h"
#include "bfstdcooc.h"

/* 
 * function to compute a sequential top-down BFS for unweighted graphs,
 * represented by sparse adjacency matrices in COOC format.   
 *  
*/
int bfs_seq_td_cooc (int *I,int *J,float *sigma,int r,int nz,int n){

    int d = 0;
    int c = 1;
    float *f;
    float *f_t;
    f =  (float *) calloc(n,sizeof(*f));
    f_t =  (float *) calloc(n,sizeof(*f_t));
    f[r] = 1.0;

    /* timing variables  */
    double initial_t;
    double delta;
    double sum_vs_t = 0;
    double spmv_t = 0;
    double assign_v_t = 0;
    double mult_vs_t = 0;
    double check_f_t = 0;
    double total_t;
    
    while (c > 0) {
      d++;
      initial_t = get_time();
      sum_vs (sigma,f,n);
      delta = get_time()-initial_t;
      sum_vs_t += delta;
      
      initial_t = get_time();
      spmv_seq_td_cooc (f,I,J,f_t,nz);
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
      c = 0;
      check_f(&c,f,n);
      delta = get_time()-initial_t;
      check_f_t += delta;
    }

    printf("\nbfs_seq_td_cooc::d = %d,r = %d \n",d,r);
    total_t =  sum_vs_t +  spmv_t + assign_v_t +  mult_vs_t + check_f_t; 
    int p_t = 1;
    if (p_t) {    
      printf("bfstdcooc_seq::f <-- f +sigma time = %lfs \n",sum_vs_t);
      printf("bfstdcooc_seq::f_t <-- A'f time = %lfs \n",spmv_t);
      printf("bfstdcooc_seq::f <-- f_t time = %lfs \n",assign_v_t);
      printf("bfstdcooc_seqc::f <-- f*(-sigma) time = %lfs \n",mult_vs_t);
      printf("bfstdcooc_seq::c <-- check (f=0)) time = %lfs \n",check_f_t);
      printf("bfstdcooc_seq::total time = %lfs \n",total_t);
    }
    
    return 0;
}//end  bfs_seq_td_cooc
