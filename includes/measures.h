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
/****************************************************************************************
*               calculate current average neighborhood vacancy (for hosts)	        *
****************************************************************************************/
double calcAverNeighborhoodVacancy(void);
/****************************************************************
*       Calculate the average investment in each cluster,       *
*       the label of the cluster with the highest               *
*       investment and its size. Return the system average      *
*       investment.                                             *
*****************************************************************/
double findHighestInvCluster(int ncl,int *labels,double *wcl,double *clwstats,ClusterFullID *clid);
/********************************************************
*  Find the number of clusters with similar investment  *
*********************************************************/
int findSimilarInvClusters(int numcl,double targetw, double eps, double *wcl);
/********************************************************
*  Calculate the difference in microbial 		*
*  composition between host @idp and its children @idk.*
*  Store result in a vector where each element		* 
*  correspond to a different event			*
********************************************************/
void storeBacDiffComp(int idp,int idk);
/**********************************************
*   store the frequency of beneficial         *
*   bacteria in a single host                 *
**********************************************/
void averInvestmentXt(void);
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
/************************************************************************
*       Stores the following cluster related measures:                  *
*       ->Highest cluster average investment (=best clluster)           *
*       ->Cluster size of the best cluster                              *
*       ->Investment of the largest cluster                             *
*       ->Size of the largest cluster                                   *
*       ->number of clusters                                            *
*       ->average cluster size and related standart deviation           *
*       ->System average investment                                     *
*       ->Number of hosts                                               *
*************************************************************************/
void clustersXt(void);
/********************************************************
*       Stores transient investment distribution among  *
*       clusters:                                       *
*       -> time                                         *
*       -> all possible average investments, <w>        *
*       ->fraction of clusters in each                  *
*         bin (int)(<w>/binsize)                        *
*********************************************************/
void clustersInvDistXt(void);
/***************************************************
*  stores the number of host events per            *
*  microbial time steps                            *
****************************************************/
void numHostEventsPerDtXt(void);
/****************************************************************
*       store the von neuman distance, fraction of pair of      *
*       hosts of the time distant at each von neumann distance, *
*        average investment and time                            *
*****************************************************************/
void rvnXtXw(void);
/****************************************************************
*       store spatial correlation for lattice ocupation         *
*       as a function of time and distance                      *
*****************************************************************/
void spatialCorrXt(void);
/****************************************************************
*       store spatial correlation for the idividual investment  *
*       as a function of time and distance                      *
*****************************************************************/
void spatialInvestmentCorrXt(void);
/********************************************************
*       Measure average difference of hosts parents     *
*       and their childrens microbial composition       *
*       over time (to compare different vertical        *
*       transmission approachs)                         *
*********************************************************/
void difMicrCompXt(void);
/************************************************************************
*  Store in @SAMPLE files the average investment over time              *
*  for multiple costs (@sample files for each cost).                    *
*  The routine changes the values of the cost, but the measuring and    *
*  storing is made by @averInvestmentXt(), which is being called in the         *
*  the evolution routine @callSysDynamics().                            *
*************************************************************************/
void averInvXtMultipleCosts(Event *event,Event *mevent);
/****************************************************************
*  Store the steady state average investment for different      *
*  costs and probability values of long range migration, plr    *
*****************************************************************/
void costXplrXw(Event *event,Event *mevent);
/****************************************************************
*  Store in @SAMPLE files the average investment                *
*  in the system as a function of the system carrying dilution  *
*****************************************************************/
void averInvXrh(Event *event,Event *mevent);
/************************************************************************
*  Store in @SAMPLE files the average investment                        *
*  in the system as a function of the system number of microbial        *
*  generations per host generations, Gh                                 *
*************************************************************************/
void averInvXgh(Event *event,Event *mevent);
/************************************************************************
*  Store in @SAMPLE files the average investment                        *
*  in the system as a function of the microbial migration rate          *
*************************************************************************/
void averInvXmb(Event *event,Event *mevent,double mbmin,double mbmax);
/************************************************************************
*  Store in @SAMPLE files the average investment                        *
*  in the system as a function of the fraction of long-range host       *
*  migration                                                            *
*************************************************************************/
void averInvXplr(Event *event,Event *mevent);
/************************************************************************
*  Store in @SAMPLE files the average investment                        *
*  in the system as a function of the cost factor payed by helpfull     *
*  bacteria                                                             *
*************************************************************************/
void averInvXcost(Event *event,Event *mevent);
/************************************************************************
*  Store in @SAMPLE files the average investment                        *
*  in the system as a function of the host migration coefficient mh     *
*************************************************************************/
void averInvXmh(Event *event,Event *mevent);
/************************************************************************
*               Store heatmap rhXmhX<w>                                 *
*       (cost X bac. migr. rate X average investment)                   *
*************************************************************************/
void costXmbXw(Event *event,Event *mevent);
/****************************************************************************************
*  		Store heatmap rhXmhX<w>asy						* 
*  	(Kh/N X host migration rate coefficient X asymptotic average investment)	*
*****************************************************************************************/
void rhXmhXw(Event *event,Event *mevent);
/************************************************************************
*  		Store heatmap mhXplrX<w>asy				* 
* 		*mh=migration rate coefficient  			*
* 		*plr=frac. of long-range host migrations 		*
* 		*<w>=asymptotic average investment			*
*************************************************************************/
void mhXplrXw(Event *event,Event *mevent);
/************************************************************************
*  Store in @SAMPLE files the average investment                        *
*  in the system as a function of the cost factor payed by helpfull     *
*  bacteria                                                             *
*************************************************************************/
void averInvXcost(Event *event,Event *mevent);
/***************************************************
*  call routines that measure and store measures   *
*  during the time loop                            *
****************************************************/
void timeMeasures(void);
#endif
