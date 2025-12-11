#ifndef RANDGEN_UFRGS_H
#define RANDGEN_UFRGS_H

#define FNORM   (2.3283064365e-10)
#define RANDOM  ( (ira[ip++]=ira[ip1++]+ira[ip2++]) ^ira[ip3++] )        
#define FRANDOM (FNORM * RANDOM)

extern unsigned zseed, ira[256];
extern unsigned char ip,ip1,ip2,ip3;

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
unsigned rand4init(void);
/*****************************************/
void Init_Random(void);
/****************************************************************
*                          Random Seed                          *
*                  Last Modified: 19/11/2012                    *
*                                                               *
* When debugging, use always the same seed, otherwise, take it  *
* as the numbers of seconds ellapsed since ?/?/1970.            *
****************************************************************/ 
unsigned long start_randomic(unsigned long semente);
/*********************************************************************
***                    Measuring Performance (time)                ***
***                     Last modified: 24/03/1999                  ***
***                                                                ***
*** Use:                                                           ***
***  cpu_time = (clock() - begin_clock) / (double) CLOCKS_PER_SEC; ***
***  real_time = difftime(time(NULL), begin_time);                 ***
*********************************************************************/ 
void start_stopwatches(void);
/**********************************************************************
***                          Create Time Table                      ***
***                    Last Modified: 24/01/2006                    ***
**********************************************************************/
void create_time_table(long *t1,long *t2,int power, int base);
/*********************************************************************
***                         Time Table 2                           ***
***                    Last Modified: 08/01/2000                   ***
***  The total number of time steps and the number of measures are ***
***  specified.                                                    ***
*********************************************************************/
void create_time_table_2(unsigned long *t1, unsigned long total_time, unsigned long measures);
/*********************************************************************
***                       Time Table (linear)                      ***
***                    Last Modified: 05/08/2011                   ***
***  The total number of time steps and the number of measures are ***
***  specified.                                                    ***
*********************************************************************/
void create_time_table_linear(unsigned long *t1, unsigned long total_time, int measures, int initial);
/*********************************************************************
***                        Time Table (log)                        ***
***                    Last Modified: 07/01/2003                   ***
***  The number of decades of times is specified.                  ***
*********************************************************************/
void create_time_table_decades(unsigned long *t1, int decades);
#endif
