#ifndef EVO_H
#define EVO_H
/****************************************************************
*     calculates the accumulated investment in host @index      *
*****************************************************************/
double calcAcumInvest(int index,SysParams *spar);
/****************************************************
* Calculates the host event rates: birth and death  *
*****************************************************/
void calcHostEvents(double *event, SysParams *spar);
/***************************************************************
*     Set microbial frequencies of the offspring of a host     *
*     located. The parent and the kid hosts are located at     *
*     sites idp and idk, respectively                          *
****************************************************************/
void setMicrKids(int idp, int idk,SysParams *spar);
/******************************************************
 *      birth of a new host                           *
 *****************************************************/
void hostBirth(int idp, int idk,SysParams *spar);
/******************************************************
*       death of a host                               *
*******************************************************/
void hostDeath(int idh);
/************************************************
*         returns gillespie's time increment    *
************************************************/
double gillespieTime(double sumprob);
/***********************************************
*         Adjust host time step                *
************************************************/
double adjustTimeStep(double *event,int size);
/**************************************************
*      Adjust host time step: simple version      *
*      Chooses the maximum vector element @maxe   *
*      and considers that the maximum probability *
*      of 2 host events is ~maxe*maxe             *
**************************************************/
double adjustTimeStepSimple(double *event,int size);
/*******************************************************
*                  host dynamics                       *
********************************************************/
void dynamicsHost(Event *event,SysParams *spar);
#endif
