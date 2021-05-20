/*
 * This source code is distributed under the terms defined on 
 * in the file bfstdcooc_main.c of this source distribution.
*/
/* 
*   Breadth first search (BFS)  
*   Single precision (float data type)
*   TurboBFS_COOC_TD:utils_bfs.h
* 
*  Functions to:
*  1) print results
*  2) check the GPU computed values of BFS shortest paths by comparing
*     with the results of the sequential BFS.
* 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 
 * function to print results
*/
int printBFS(int *IC,int *JC,float *sigma,float *sigma_hgpu,int nz, int n){	      

    int i,m1 = 100;

    /*printf("CoocA.IC  CoocA.JC \n");    
    if (n > m1){
       for (i=0; i<m1; i++){
	  printf("%d %d\n", i, IC[i], JC[i]);
	}
	for (i=nz-m1; i<nz; i++){
	  printf("%d %d\n", i, IC[i], JC[i]);
	}
    }else{
	for (i=0; i<nz; i++){
	  printf("%d %d\n", i, IC[i], JC[i]);
	}
    }*/

    printf("\nsigma   sigma_hgpu \n");
    if (n > m1){
	for (i=0; i<m1; i++){
	  printf("%d,%lf,%lf\n", i,sigma[i],sigma_hgpu[i]);
	}
	for (i=n-m1; i<n; i++){
	  printf("%d,%lf,%lf\n", i,sigma[i],sigma_hgpu[i]);
	}
    }else{
	for (i=0; i<n; i++){
	  printf("%d,%lf,%lf\n", i,sigma[i],sigma_hgpu[i]);
	}
    }
    return 0;
}// end printBFS

/******************************************************************************/
/* 
 * function to check the GPU computed values of BFS shortest paths by comparing
 * with the results of the sequential BFS. 
*/
int bfs_check(float *sigma,float *sigma_hgpu,int n){

    int i;
    int check = 1;
    float epsilon = 1.0e-5;
    int count = 0;
    
    for (i=0; i<n; i++){
      if (sigma[i] > n || sigma_hgpu[i] > n  ){
	  printf("error: value on computed sigma vector greater than the number of vertices \n");
	  printf("i = %d, sigma[%d] = %lf, sigma_hgpu[%d] = %lf \n", i, i, sigma[i], i, sigma_hgpu[i]);
 	  return 1;
      }
      if (fabs (sigma[i]-sigma_hgpu[i]) > epsilon){
	check = 0;
        count++;      
      }
    }
    if(check){
      printf("values of the GPU computed sigma vector are correct \n");
    }else{
      printf("error in %d values of the GPU computed sigma_hgpu vector\n",count);
    }	
  
    return 0;
}//end bfs_check


