/*
 * This source code is distributed under the terms defined on 
 * in the file bfstdcsc_main.c of this source distribution.
*/
/* 
*   Breadth first search (BFS)  
*   Single precision (float data type)
*   TurboBFS_CSC_TD:utils_bfs.c
* 
*  Functions to:
*  1) print results
*  2) check the GPU computed values of BFS shortest paths by comparing
*    with the results of the sequential BFS.
*  3) check the GPU computed values of the S vector by comparing
*     with the results of the sequential computed values of S.
* 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 
 * function to print results
*/
int printBFS(int *IC,int *CP,int *S,int *S_hgpu,float *sigma,float *sigma_hgpu,
	      int nz, int n){

    int i,m1 = 100;

    /*printf("CscA.IC\n");    
    if (n > m1){
       for (i=0; i<m1; i++){
	  printf("%d %d\n", i, IC[i]);
	}
	for (i=nz-m1; i<nz; i++){
	  printf("%d %d\n", i, IC[i]);
	}
    }else{
	for (i=0; i<nz; i++){
	  printf("%d %d\n", i, IC[i]);
	}
    }
    printf("CscA.CP \n");
    if (n > m1){
	for (i=0; i<m1; i++){
	  printf("%d %d\n", i, CP[i]);
	}
	for (i=n-m1; i<n+1; i++){
	  printf("%d %d\n", i, CP[i]);
	}
    }else{
	for (i=0; i<n+1; i++){
	  printf("%d %d\n", i, CP[i]);
	}
	}*/

    printf("\nS   S_hgpu \n");
    if (n > m1){
	for (i=0; i<m1; i++){
	  printf("%d,%d,%d\n", i,S[i],S_hgpu[i]);
	}
	for (i=n-m1; i<n; i++){
	  printf("%d,%d,%d\n", i,S[i],S_hgpu[i]);
	}
    }else{
	for (i=0; i<n; i++){
	  printf("%d,%d,%d\n", i,S[i],S_hgpu[i]);
	}
    }
    
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

/******************************************************************************/
/* 
 * function to check the GPU computed values of the S vector by comparing
 * with the results of the sequential computed values of S. 
*/
int S_check(int *S,int *S_hgpu,int n){

    int i;
    int check = 1;
    int epsilon = 1.0e-5;
    int count = 0;
    
    for (i=0; i<n; i++){
      if (fabs (S[i]-S_hgpu[i]) > epsilon){
	check = 0;
        count++;      
      }
    }
    if(check){
      printf("values of the GPU computed S vector are correct \n");
    }else{
      printf("error in %d values of the GPU computed S_hgpu vector\n",count);
    }	
  
    return 0;
}//end S_check
