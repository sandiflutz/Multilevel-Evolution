/*randgen_ufrgs.c*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include"randgen_ufrgs.h"

static clock_t   begin_cpu_time;                 /* Time variables */
static time_t    begin_real_time;
unsigned long randomize, seed;                    /* Random number */
unsigned zseed, ira[256];
unsigned char ip,ip1,ip2,ip3;
/********************************************************************
*            Random Number Generator by Parisi & Rapuano            *
*                  Last Modified: 18/05/2000                        *
*                                                                   *
*  First, the function start_randomic() should be called to create  *
*  the seed (odd number):                                           *
*                 seed = start_randomic();                          *
*                                                                   *
*  Return: start_randomic -> unsigned long int                      *
*          FRANDOM -> double in [0,1)                               *
********************************************************************/
unsigned rand4init(void){
     unsigned long long y;
   
     y = (zseed*16807LL);
     zseed = (y&0x7fffffff) + (y>>31);
     if (zseed&0x80000000){
	     zseed = (zseed&0x7fffffff) + 1;
     }

     return zseed;
}
/********************************************************/
void Init_Random(void){
     int i;
   
     ip=128;
     ip1=ip-24;
     ip2=ip-55;
     ip3=ip-61;
   
     for (i=ip3; i<ip; i++){
	     ira[i] = rand4init();
     }
}
/****************************************************************
*                          Random Seed                          *
*                  Last Modified: 19/11/2012                    *
*                                                               *
* When debugging, use always the same seed, otherwise, take it  *
* as the numbers of seconds ellapsed since ?/?/1970.            *
****************************************************************/ 
unsigned long start_randomic(unsigned long semente){
	
	if (semente==0) semente = (unsigned long) time(NULL);      /* random seed */
	if (semente%2==0) ++semente; /* odd number */
	#ifdef DEBUG
	 semente = 1016476867;
	#endif
	zseed = semente;
	Init_Random();
	return semente;
}
/*********************************************************************
***                    Measuring Performance (time)                ***
***                     Last modified: 24/03/1999                  ***
***                                                                ***
*** Use:                                                           ***
***  cpu_time = (clock() - begin_clock) / (double) CLOCKS_PER_SEC; ***
***  real_time = difftime(time(NULL), begin_time);                 ***
*********************************************************************/ 
void start_stopwatches(void){
	begin_cpu_time = clock();
	begin_real_time = time(NULL);
	return;
}
/**********************************************************************
***                          Create Time Table                      ***
***                    Last Modified: 24/01/2006                    ***
**********************************************************************/
void create_time_table(long *t1,long *t2,int power, int base){
	long i,j;
   
	t1[0] = 1;
	for (i=1; i<power; i++){
		t1[i] = t1[i-1]*base;
	}

	if (t2 != NULL){
    		for (i=0; i<power; i++){
        		for (j=0; j<power; j++){
            			*(t2+i*power+j) = t1[i] + t1[j];
				if (*(t2+i*power+j) >= base*t1[power-1]) *(t2+i*power+j) =-1;
   
			}
		}
	}

	return;
}                    
/*********************************************************************
***                         Time Table 2                           ***
***                    Last Modified: 08/01/2000                   ***
***  The total number of time steps and the number of measures are ***
***  specified.                                                    ***
*********************************************************************/
void create_time_table_2(unsigned long *t1, unsigned long total_time, unsigned long measures){
	unsigned long i,k;
	double temp;
   
	t1[0] = 0;
	temp = pow((double) total_time,1.0/(measures-1));
	k=0;
	for (i=1; i<measures; ++i){
	   	t1[i] = (int) pow(temp, (double) i);
	   	if (t1[i]<=k) t1[i]=k+1;
	   	k = t1[i];
    	}

	return;
}
/*********************************************************************
***                       Time Table (linear)                      ***
***                    Last Modified: 05/08/2011                   ***
***  The total number of time steps and the number of measures are ***
***  specified.                                                    ***
*********************************************************************/
void create_time_table_linear(unsigned long *t1, unsigned long total_time, int measures, int initial){
	unsigned long i;
	double temp;
   
	t1[0] = 0;
	temp = ((double) total_time)/measures;
	if (initial==0){
		t1[1] = temp;
	}else{ 
		t1[1] = initial; /* thermalization time */
	}
	for (i=2; i<=measures; ++i){
     		t1[i] = t1[i-1] + (int) temp;
	}

	return;
}
/*********************************************************************
***                        Time Table (log)                        ***
***                    Last Modified: 07/01/2003                   ***
***  The number of decades of times is specified.                  ***
*********************************************************************/
void create_time_table_decades(unsigned long *t1, int decades){
	unsigned long i=1,j,k=1;

	t1[0] = 0;
	while (i <= decades){
		t1[k] = (int) pow(10,i-1);
		++k;
		for (j=2; j<10; ++j){
			t1[k] = t1[k-1] + (int) pow(10,i-1);
			++k;
		}
		++i;
	}
	t1[k] = (int) pow(10,i-1);

	return;
}
