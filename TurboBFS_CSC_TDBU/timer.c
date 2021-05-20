
/* 
* 
*  This program provides a timer. 
* 
*/

#include<time.h>
#include "timer.h"
#include <sys/time.h>


struct timeval start;

void reset_timer(){
  gettimeofday(&start, NULL);
}

double get_time(){
  struct timeval end;
  gettimeofday(&end, NULL);
  return ((double)(end.tv_sec-start.tv_sec) + (double)(end.tv_usec/1000000.0 - start.tv_usec/1000000.0));
} 


