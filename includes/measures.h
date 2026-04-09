#ifndef MEASURES_H
#define MEASURES_H
/****************************************************************
*     allocate memory for global arrays and structs        	*
*     related to time measuraments. Initialize time 		*
*     measurement parameters					*
****************************************************************/
void allocateMemTM(void);
/************************************************
*       Free Allocated Memory for time          *
*       measurements                            *
*************************************************/
void freeMemTM(void);
/******************************************
*       Open Global Files                 *
*******************************************/
void openFiles(void);
/**********************************************
*               Close Global Files            *
***********************************************/
void closeFiles(void);
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
double calcAverInv(void);
/********************************************************
*  Calculate the difference in microbial 		*
*  composition between host @idp and its children @idk.*
*  Store result in a vector where each element		* 
*  correspond to a different event			*
********************************************************/
void storeBacDiffComp(int idp,int idk);
/************************************************
*   store the current average investment level  *
*   in an isolate host                          *
**************************************************/
void densB1Xt(void);
/********************************************************
* Stores the mean fraction of each type of bacteria in  *
* the system                                            *
*********************************************************/
void meanFracXt(void);
/**********************************************
*   store the frequency of beneficial         *
*   bacteria in a single host                 *
**********************************************/
void averInvestmentXt(void);
/********************************************************
*       Store measures related to vancancy              *
*       frequency as functions  of time                 *
*       -> average group vacancy freq.,<rho_e>          *
*       -> standart deviation of rho_e                  *
*       -> (Kh-H(t))/N=Rho_e-Rho_ekh                    *
*       (where Rho_e is the system vacancy              *
*       freq. and Rho_ekh is the vacancy freq.          *
*       related to the carrying capacity=(N-Kh)/N       *
*********************************************************/
void emptyFreqXt(void);
/**************************************
*  snapshot of host network           *
*  The colors indicate the acumulated *
*  investment of each host            *
***************************************/
void save_config(void);
/***********************************************
 * store investment distribution among and     *
 * create a gnuplot script to create graphics  *
 ***********************************************/
void invDistXt(void);
/***************************************************
*  stores the number of host events per            *
*  microbial time steps                            *
****************************************************/
void numHostEventsPerDtXt(void);
/***************************************************
*  calculate spatial corretation: both horizontal  *
*  and vertical for a specific distance,           *
*  considering                                     *
****************************************************/
void spatialCorrXt(void);
/********************************************************
*       Measure average difference of hosts parents     *
*       and their childrens microbial composition       *
*       over time (to compare different vertical        *
*       transmission approachs)                         *
*********************************************************/
void difMicrCompXt(void);
/***************************************************
*  storing average host generation time in a file  *
****************************************************/
void genHostTime(void);
/****************************************************************
*  Store in @SAMPLE files the average investment                *
*  in the system as a function of the system carrying dilution  *
*****************************************************************/
void averInvXrh(Event *event,Event *mevent);
/***************************************************
*  call routines that measure and store measures   *
*  during the time loop                            *
****************************************************/
void timeMeasures(void);
#endif
