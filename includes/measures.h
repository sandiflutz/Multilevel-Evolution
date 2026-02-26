#ifndef MEASURES_H
#define MEASURES_H
/******************************************
*       Open Global Files                 *
*******************************************/
void openFiles(TimeMeasures *meas);
/**********************************************
*               Close Global Files            *
***********************************************/
void closeFiles(TimeMeasures *meas);
/**********************************************
 * calculate investment density per Host     *
 **********************************************/
void calcInvDens(double *densInvH);
/**********************************************
 * calculate investment frequency per Host    *
 **********************************************/
void calcInvFreq(double *freqInvH);
/**********************************************************************
 * calculate investment distribution among                            *
 * hosts                                                              *
 **********************************************************************/
void calcInvDist(double binsize, double *hist_inv, double *freqInvH); 
/************************************************
*   store the current average investment level  *
*   in an isolate host                          *
**************************************************/
void densB1Xt(TimeMeasures *meas);
/********************************************************
* Stores the mean fraction of each type of bacteria in  *
* the system                                            *
*********************************************************/
void meanFracXt(TimeMeasures *meas);
/**********************************************
*   store the frequency of beneficial         *
*   bacteria in a single host                 *
**********************************************/
void averInvestmentXt(TimeMeasures *meas);
/**************************************
*  snapshot of host network           *
*  The colors indicate the acumulated *
*  investment of each host            *
***************************************/
void save_config(TimeMeasures *meas);
/***********************************************
 * store investment distribution among and     *
 * create a gnuplot script to create graphics  *
 ***********************************************/
void invDistXt(TimeMeasures *meas);
/***************************************************
*  stores the number of host events per            *
*  microbial time steps                            *
****************************************************/
void numHostEventsPerDtXt(TimeMeasures *meas);
/***************************************************
*  call routines that measure and store measures   *
*  during the time loop                            *
****************************************************/
void measures(TimeMeasures *meas);
#endif
