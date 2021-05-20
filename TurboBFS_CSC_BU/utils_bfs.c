/*
 * This source code is distributed under the terms defined on 
 * in the file bfsbucsc_main.c of this source distribution.
*/
/* 
*   Breadth first search (BFS)  
*   Single precision (float data type)
*   TurboBFS_CSC_BU:utils_bfs.c
* 
*  This program:
*  1) prints results
*  2) check the GPU computed values of the S vector by comparing
*     with the results of the sequential computed values of S.
* 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 
 * function to print results
*/
int printBFS(int *I,int *CP,int *S,int *S_hgpu,float *sigma,int nz, int n){
	      
    int i,m1 = 100;

    /*printf("CscA.IC\n");    
    if (n > m1){
       for (i=0; i<m1; i++){
	  printf("%d %d\n", i, I[i]);
	}
	for (i=nz-m1; i<nz; i++){
	  printf("%d %d\n", i, I[i]);
	}
    }else{
	for (i=0; i<nz; i++){
	  printf("%d %d\n", i, I[i]);
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
    
    printf("\nsigma_seq  \n");
    if (n > m1){
	for (i=0; i<m1; i++){
	  printf("%d,%lf\n", i,sigma[i]);
	}
	for (i=n-m1; i<n; i++){
	  printf("%d,%lf\n", i,sigma[i]);
	}
    }else{
	for (i=0; i<n; i++){
	  printf("%d,%lf\n", i,sigma[i]);
	}
    }
    return 0;
}// end printBFS

/******************************************************************************/
/* 
 * function to check the GPU computed values of the S vector by comparing
 * with the results of the sequential computed values of S. 
*/
int S_check(int *S,int *S_hgpu,int n){

    int i;
    int check = 1;
    int epsilon = 1.0e-3;
    int count = 0;
    
    for (i=0; i<n; i++){
      if (fabs(S[i]-S_hgpu[i]) > epsilon){
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
