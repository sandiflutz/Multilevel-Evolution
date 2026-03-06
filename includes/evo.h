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
/****************************************************************
 *      Returns the number of host time steps per bacteria      *
 *      time step                                               *
 ****************************************************************/
int hostNTSPerBacNTS(Event *event);
/****************************************************************
*       Host Layer Evolution (for a time interval=Dt_ref) for   *
*       the complete graph version                              *
*****************************************************************/
int evolveHostCG(int dnumsteps,Event *event, TimeMeasures *meas);
/****************************************************************
*       Host Layer Evolution (for a time interval=Dt_ref) for   *
*       the lattice version                                     *
*****************************************************************/
void evolveHostSL(int dnumsteps,Event *event, TimeMeasures *meas);
#endif
