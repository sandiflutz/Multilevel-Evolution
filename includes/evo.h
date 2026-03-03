#ifndef EVO_H
#define EVO_H
/****************************************************************
*     calculates the accumulated investment in host @index      *
*****************************************************************/
double calcAcumInvest(int index);
/****************************************************
* Calculates the host event rates: birth and death  *
*****************************************************/
void calcHostEvents(Event *event);
/***************************************************************
*     Set microbial frequencies for the offspring of hosr @idp *
*     Bacteria types and their frequencies are randomly        *
*     selected using a normal distrution around the bacteria   *
*     frequencies on the parent                                *
****************************************************************/
void setMicrKidsNorm(int idp, int idk);
/***************************************************************
*     Set microbial frequencies for the offspring of hosr @idp *
*     Bacteria types and their frequencies are randomly        *
*     selected using a poisson distribution for the number of  *
*     bacteria types in the offspring                          *
****************************************************************/
void setMicrKidsPoiss(int idp, int idk);
/******************************************************
 *      birth of a new host                           *
 *****************************************************/
void hostBirth(int idp, int idk);
/******************************************************
*       death of a host                               *
*******************************************************/
void hostDeath(int idh);
/****************************************
*       moviment of hosts:              *
*       2 neighbors exchange places     *
*****************************************/
void hostMoviment(int id1,int id2);
/************************************************
*         returns gillespie's time increment    *
************************************************/
double gillespieTime(double sumprob);
/**************************************************
*      Adjust host time step: the probability of  *
*      2 host events in a host timestep os <0.01  *
***************************************************/
double adjustTimeStep(double maxprob);
/*******************************************************
*                  host dynamics                       *
********************************************************/
void dynamicsHost(int type_event,int idh,DynList *lhost, int *ilhost,TimeMeasures *meas);
/*************************************************************
*    Host Layer Evolution (for a time interval=Dt_ref):      *
*    This version uses a dynamical value for the time step   *
*    of the host layer of dt<=Dt_ref                         *
*    (original paper version)                                *
**************************************************************/
int evolveHostDtH(Event *event,DynList *listh_tmp, int *inverselisth_tmp,TimeMeasures *meas);
/*************************************************************
*    Host Layer Evolution (for a time interval=Dt_ref):      *
*    This version uses a tau-leaping method                  *
**************************************************************/
void evolveHostTLP(Event *event,DynList *listh_tmp, int *inverselisth_tmp,TimeMeasures *meas);
/*************************************************************
*    Host Layer Evolution (for a time interval=dtE):         *
*    This version uses mont carlo steps                      *
**************************************************************/
void evolveHostMCS(Event *event,DynList *listh_tmp, int *inverselisth_tmp,TimeMeasures *meas);
#endif
