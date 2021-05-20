
/*
 *   Functions to:
 *   1) compute the parameters(maximum, mean, and standard
 *      deviation) of the degree, out degree and in degree
 *      vectors
 *   2) compute the degree vector for undirected graphs.
 *
 *
 */

#ifndef DEGREE_H
#define DEGREE_H


/* 
 * function to compute the degree vector of an undirected graph represented by a 
 * symmetric sparse adjacency matrix.
 */
int degree_ug (int *out_degree,int *in_degree,int *degree,int N);

/******************************************************************************/
/* 
 * function to compute the value and index of the maximum element of the degree
 * vector.
 *
*/   
int degree_max (int *degree,int *maxDegree,int *idxMaxDegree,int N);

/******************************************************************************/
/*  function to compute: the mean value, the variance, and the standard deviation 
 *  of the values in the degree vector.  
 */    
int degree_sta (int *degree,double *meanDegree,double *sdDegree,int N);

/******************************************************************************/
/* 
 * function to compute the degree distribution vector of a graph represented by a 
 * sparse adjacency matrix.
 */
int degree_dist (int *degree,int *degreeDist,int N);


/******************************************************************************/
/* 
 * function to compute the relative degree distribution vector of a graph 
 * represented by a sparse adjacency matrix.
 */
int degree_dist_rel (int *degreeDist,double *degreeDistRel,int maxDegree,int N);

/******************************************************************************/
/* 
 * function to compute the accumulated degree distribution and relative degree 
 * distribution vectors of a graph represented by a sparse adjacency matrix.
 * 
 */
int degree_dist_ac (int *degreeDist,double *degreeDistRel,int *degreeDistAc,
		    double *degreeDistRelAc,int maxDegree);


#endif
