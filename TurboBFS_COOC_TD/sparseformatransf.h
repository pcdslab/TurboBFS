/*
 * This source code is distributed under the terms defined on 
 * in the file bfstdcooc_main.c of this source distribution.
 */
/* 
 *  Breadth first search (BFS)
 *  Single precision (float data type)
 *  TurboBFS_COOC_TD:sparseformatransf.h
 *   
 *  Functions:
 *  1) to transform asymmetric sparse matrices in Compressed Sparse Column (CSC) 
 *     format to Compressed Sparse Row (CSR) format,
 *  2) symmetric sparse matrices on Symmetric Sparse Column (SSC) 
 *     format to coordinate ordered by column (COOC) format, 
 *  3) to compute CP and RP arrays.
 *
 */

#ifndef SPARSE_FORMATRANSF_H
#define SPARSE_FORMATRANSF_H


////////////////////////////////////////////////////////////////////////////////
/* 
 * function to transform sparse matrices, representing unweighted graphs, in 
 * Compressed Sparse Column (CSC) format to Compressed Sparse Row (CSR) 
 * format.
 * input:  IC,CP: CSC format, RP: CSR format
 * output: JR: CSR format
 */
int CSCCSR_uw (int *IC,int *CP,int *RP,int *JR,int n);
////////////////////////////////////////////////////////////////////////////////
/* 
 * function to transform symmetric sparse matrices, representing undirected, 
 * unweighted graphs, in Symmetric Sparse Column (SSC) format to Coordinated 
 * Ordered by Column (COOC) format.
 *
 *  input:  ICLT,CPLT,RPLT:SSC format, RP:CSR format, CP:CSC format
 *  output: ICOOC,JCOOC:COOC format 
 */
int SSCCOOC_uw (int *ICLT,int *CPLT,int *RPLT,int *RP,int *CP,int *ICOOC,
		int *JCOOC,int n);

/******************************************************************************/
/* 
 * function to compute the CPLT, RPLT, CP and RP arrays (size = N+1) for the 
 * SSC, CSC and CSR sparse formats for sparse adjacency matrices representing 
 * undirected graphs.
 */
int CP_RP_ug (int *out_degree,int *in_degree,int *degree,int *CPLT,int *RPLT,
	      int *CP,int *RP,int n);

/******************************************************************************/
/* 
 * function to compute the CP and RP arrays (size = N+1) for the CSC and CSR 
 * sparse formats representing sparse adjacency matrices of directed graphs
 */
int CP_RP_dg (int *out_degree,int *in_degree,int *CP,int *RP,int n);

#endif
