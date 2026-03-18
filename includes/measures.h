#ifndef MEASURES_H
#define MEASURES_H
/****************************************************************
*     allocate memory for global arrays and structs        	*
*     related to time measuraments. Initialize time 		*
*     measurement parameters					*
****************************************************************/
void allocateMemTM(TimeMeasures *meas);
/************************************************
*       Free Allocated Memory for time          *
*       measurements                            *
*************************************************/
void freeMemTM(TimeMeasures *meas);
/******************************************
*       Open Global Files                 *
*******************************************/
void openFiles(TimeMeasures *meas);
/**********************************************
*               Close Global Files            *
***********************************************/
void closeFiles(TimeMeasures *meas);
/****************************************************************
 *      calculate investment density per Host:                  *
 *      densInvH[host]=sum_type(bac[host][type]*inv[type])      *
 ****************************************************************/
void calcInvDens(double *densInvH);
/************************************************************************
*      calculate investment frequency per Host:                        	*
*      freqInvH[host]=sum_type(bac[host][type]*inv[type])/micr[host]   	*
************************************************************************/
void calcInvFreq(double *freqInvH);
/****************************************************************
*      calculate investment distribution among hosts           *
****************************************************************/
void calcInvDist(double binsize,double *hist_inv,double *freqInvH);
/****************************************************************************************
*       calculate current average investment in the system:                             *
*       averInv=sum_host(sum_type(bac[host][type]))/total amount of bac. in the system  *
****************************************************************************************/
double calcAverInv(double *densInvH);
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
*  calculate spatial corretation: both horizontal  *
*  and vertical for a specific distance,           *
*  considering                                     *
****************************************************/
void spatialCorrXt(TimeMeasures *meas);
/********************************************************
 *  store average difference in microbial composition   *
 *  between host @idp and its children @idk             *
 ********************************************************/
void storeBacDiffComp(int idp,int idk);
/********************************************************
*       Measure average difference of hosts parents     *
*       and their childrens microbial composition       *
*       over time (to compare different vertical        *
*       transmission approachs)                         *
*********************************************************/
void difMicrCompXt(TimeMeasures *meas);
/***************************************************
*  calculating average host generation time        *
****************************************************/
void genHostTime(TimeMeasures *meas);
/***************************************************
*  call routines that measure and store measures   *
*  during the time loop                            *
****************************************************/
void measures(TimeMeasures *meas);
#endif
