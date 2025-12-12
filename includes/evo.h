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
/**************************************************
*      Adjust host time step: the probability of  *
*      2 host events in a host timestep os <0.01  *
***************************************************/
double adjustTimeStep(double maxprob);
/*******************************************************
*                  host dynamics                       *
********************************************************/
void dynamicsHost(Event *event,SysParams *spar);
#endif
