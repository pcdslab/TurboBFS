/*
 * This source code is distributed under the terms defined on 
 * in the file bfstdcooc_main.c of this source distribution.
 */
/* 
 *  Breadth first search (BFS)
 *  Single precision (float data type)
 *  TurboBFS_COOC_TD:sparseformatransf.c
 *   
 *  Functions:
 *  1) to transform asymmetric sparse matrices in Compressed Sparse 
 *     Column (CSC) format to Compressed Sparse Row (CSR) format,
 *  2) symmetric sparse matrices on Symmetric Sparse Column (SSC) 
 *     format to coordinate ordered by column (COOC) format, 
 *  3) to compute CP and RP arrays.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "sparseformatransf.h"


/* 
 * function to transform sparse matrices, representing unweighted graphs, in 
 * Compressed Sparse Column (CSC) format to Compressed Sparse Row (CSR) 
 * format.
 * input:  IC,CP: CSC format, RP: CSR format
 * output: JR: CSR format
 */
int CSCCSR_uw (int *IC,int *CP,int *RP,int *JR,int n){

  int i,j,k,next;

  //compute JR
  for (j=0; j<n; j++){
    for (k=CP[j]; k<CP[j+1]; k++){
      i = IC[k];
      next = RP[i];
      JR[next] = j;
      RP[i] = next +1;
    }
  }

  //reshift RP
  for (i=n-1; i>=0; i--){
    RP[i+1] = RP[i];
  }
  RP[0] = 0;

  return 0;
}//end CSCCSR_uw

/******************************************************************************/
/* 
 * function to transform symmetric sparse matrices, representing undirected, 
 * unweighted graphs, in Symmetric Sparse Column (SSC) format to Coordinated 
 * Ordered by Column (COOC) format.
 *
 *  input:  ICLT,CPLT,RPLT:SSC format, RP:CSR format, CP:CSC format
 *  output: ICOOC,JCOOC:COOC format 
 */
int SSCCOOC_uw (int *ICLT,int *CPLT,int *RPLT,int *RP,int *CP,int *ICOOC,
		int *JCOOC,int N){

  int i,k,kr;

  //compute ICOOC for lower triangular matrix
  CSCCSR_uw (ICLT,CPLT,RPLT,ICOOC,N);

  //transform SSC format to COOC format

  for (i=N-1; i>=0; i--){
    kr = RP[i];
    for (k=RPLT[i]; k<RPLT[i+1]; k++){
      ICOOC[kr] = ICOOC[k];
      kr = kr + 1;
    }
  }

  for (i=0; i<N; i++){
    for (k=CPLT[i]; k<CPLT[i+1]; k++){
      kr = k + RPLT[i+1];
      ICOOC[kr] = ICLT[k];
    }
  }

  //compute JCOOC array from CP array
  for (i=0; i<N; i++){
    for (k=CP[i]; k<CP[i+1]; k++){
      JCOOC[k] = i;
    }
  }

  return 0;
}// end SSCCOOC_uw

/******************************************************************************/
/* 
 * function to compute the CPLT, RPLT, CP and RP arrays (size = N+1) for the 
 * SSC, CSC and CSR sparse formats for sparse adjacency matrices representing 
 * undirected graphs.
 */
int CP_RP_ug (int *out_degree,int *in_degree,int *degree,int *CPLT,int *RPLT,
	      int *CP,int *RP,int n) {
  int i;

  for (i=0; i<=n; i++){
    CPLT[i+1] = CPLT[i] + in_degree[i];
    RPLT[i+1] = RPLT[i] + out_degree[i];
    CP[i+1] = CP[i] + degree[i];
    RP[i+1] = CP[i+1];
  }

  return 0;
}// end CP_RP_ug

/******************************************************************************/
/* 
 * function to compute the CP and RP arrays (size = N+1) for the CSC and CSR 
 * sparse formats representing sparse adjacency matrices of directed graphs
 */
int CP_RP_dg (int *out_degree,int *in_degree,int *CP,int *RP,int n){
  int i;

  for (i=0; i<=n; i++){
    CP[i+1] = CP[i] + in_degree[i];
    RP[i+1] = RP[i] + out_degree[i];
  }

  return 0;
}// end CP_RP_dg





