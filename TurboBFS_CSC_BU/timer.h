
/* 
* 
*  This program provides a timer. 
* 
*/

#ifndef TIMER_H_
#define TIMER_H_

#define TIMER_TYPE CLOCK_PROCESS_CPUTIME_ID

void reset_timer();
double get_time();

#endif /* TIMER_H_ */
