/*
 * This source code is distributed under the terms defined on 
 * in the file bfstdcsc_main.c of this source distribution.
 */
/* 
 *  Breadth first search (BFS)
 *  Single precision (float data type)
 *  TurboBFS_CSC_TD:sparseformatransf.h
 *
 *  Functions to:
 *  1)  transform asymmetric sparse matrices in Compressed Sparse Column (CSC) 
 *      format to Compressed Sparse Row (CSR) format,
 *   2) transform symmetric sparse matrices on Symmetric Sparse Column (SSC) 
 *      format to Compressed Sparse  Row (CSC) format,
 *   3) compute CP and RP arrays.
 *
 */
#include <stdio.h>
#include <stdlib.h>

/* 
 * function to transform sparse matrices, representing unweighted graphs,in 
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
 * unweighted graphs, in Symmetric Sparse Column (SSC) format to Compressed Sparse 
 * Column (CSC) format.
 *
 *  input:  ICLT,CPLT,RPLT:SSC format, RP:CSR format, CP:CSC format
 *  output: IC:CSC format
 */
int SSCCSC_uw (int *ICLT,int *CPLT,int *RPLT,int *RP,int *CP,int *IC,int n){

  int i,k,kr;

  //compute IC
  CSCCSR_uw (ICLT,CPLT,RPLT,IC,n);

  //transform SSC format to CSC format
  for (i=n-1; i>=0; i--){
    kr = RP[i];
    for (k=RPLT[i]; k<RPLT[i+1]; k++){
      IC[kr] = IC[k];
      kr = kr + 1;
    }
  }

  for (i=0; i<n; i++){
    for (k=CPLT[i]; k<CPLT[i+1]; k++){
      kr = k + RPLT[i+1];
      IC[kr] = ICLT[k];
    }
  }

  return 0;
}// end SSCCSC_uw

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





