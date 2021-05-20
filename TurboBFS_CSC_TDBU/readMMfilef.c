
/*
 *   Read Matrix Market program based on the ANSI C library for Matrix Market I/O
 *   Single precision (float data type)
 *  
 *   This program reads a Matrix Market format file for sparse adjacency matrices 
 *   representing weighted or unweighted graphs,and computes the in-degree and 
 *   out-degree vectors.   
 *
 *   NOTES:
 *
 *   1) Matrix Market files are always 1-based, i.e. the index of the first
 *      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
 *      OFFSETS ACCORDINGLY offsets accordingly when reading and writing 
 *      to files.
 *
 *   2) ANSI C requires one to use the "l" format modifier when reading
 *      double precision floating point numbers in scanf() and
 *      its variants.  For example, use "%lf", "%lg", or "%le"
 *      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "readMMfilef.h"

/* 
 * function to read  a Matrix Market file for a weighted graph and to 
 * compute the in-degree and out-degree vectors.
 * 
*/
int readMMfile_w (FILE *f,int *IC,int *JC,float *valC,int *out_degree,
		  int *in_degree,int nz){
 
    int i;

    for (i=0; i<nz; i++){
	if(fscanf(f, "%d %d %g\n", &IC[i], &JC[i], &valC[i]) == 3){
	   IC[i]--;  /* adjust from 1-based to 0-based */
           JC[i]--;
	}            
      out_degree[IC[i]]++;
      in_degree[JC[i]]++;  
    }
   
    return 0;
}// end readMMfile_w

/******************************************************************************/
/* 
 * function to read  a Matrix Market file for a unweighted graph and to 
 * compute the in-degree and out-degree vectors.
 * 
*/
int readMMfile_uw (FILE *f,int *IC,int *JC,int *out_degree,int *in_degree,int nz){

    int i;

    for (i=0; i<nz; i++){
	if(fscanf(f, "%d %d\n", &IC[i], &JC[i]) == 2){
	   IC[i]--;  /* adjust from 1-based to 0-based */
           JC[i]--;
	}            
      out_degree[IC[i]]++;
      in_degree[JC[i]]++;  
    }
   
    return 0;
}// end readMMfile_uw



