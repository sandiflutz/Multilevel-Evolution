#ifndef EVO_H
#define EVO_H
/****************************************************************
*     calculates the accumulated investment in host @index      *
*****************************************************************/
double calcAcumInvest(int index);
/********************************************************
*	Set host birth-death rates  			*
*********************************************************/
void setIndividualHostRates(Event *event);
/****************************************************************
*   For non well-mixed cases:                                   *
*   each host event rate, associated to host @i, is             *
*   substituted by the average host rate in the group centered  *
*   on @i                                                       *
*****************************************************************/
void setGrRates(Event *event,Event *mevent);
/****************************************************************
*       set cumulative rates for host events        		*
*****************************************************************/
void setCumulativeRates(Event *event);
/****************************************************************
*     Set microbial frequencies for the offspring of host @idp 	*
*     Bacteria types and their frequencies are randomly        	*
*     selected using a normal distribution around the bacteria	*
*     frequencies on the parent.                                *
****************************************************************/
void setMicrKidsNorm(int idp, int idk);
/****************************************************************
*     Set microbial frequencies for the offspring of host @idp  *
*     Bacteria types and their frequencies are randomly         *
*     selected using a normal distribution around the bacteria	*
*     frequencies on the parent.                                *
*     -In this version, the amount of                           *
*     bacteria passed to the children are disappear from        *
*     the parent.                                               *
*     -The amount of bacteria passed is a fixed frequency Fp    *
*     of the amount of bacteria in the parent.                  *
*     -Types of bacteria, in the parent, that have a frequency  *
*     of fp<sigma are not included the kids microbiome          *
*****************************************************************/
void setMicrKidsNormPass(int idp, int idk);
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
*       Update group vacancy frequency centered on site @idh    *
*       due to a change in its state:                           *
*       if site @idh is occupied, means there was a decrease    *
*       in the number of empty sites in the group, otherwise,   *
*       there was an increase.                                  *
*****************************************************************/
void updateEmptySpaceGrFreq(int idh);
/********************************************************
*       Choose site for migration:                      *
*       states of sites @idm and one chosen are         *
*       exchanged. Choice randomly chooses a            *
*       neighboring site, according to probabilities    *
*       that depend on vacancy (VIZ=#of neighbors):     *
*               -prob[i]=1/VIZ, if neighbor is empty,   *
*               -prob[i]=rho_e[i]/VIZ, otherwise        *
*********************************************************/
int chooseMigSite(int idm);
/****************************************************************
*       Kill hosts that have a microbiome extremely low         *
*****************************************************************/
void killHostWithoutMicr(void);
/************************************************
*       Dynamics for host migration events      *
*       that happen in a Dt_ref (=microbial     *
*       time step)                              *
*************************************************/
void hostMigrationDynamics(int dnumsteps,Event *mevent);
/****************************************************************
*       Host Layer Evolution (for a time interval=Dt_ref) for   *
*       the complete graph version                              *
*****************************************************************/
void evolveHostCG(int dnumsteps,Event *event);
/****************************************************************
*       Host Layer Evolution (for a time interval=Dt_ref) for   *
*       the lattice version                                     *
*****************************************************************/
void evolveHostSL(int dnumsteps,Event *event);
/****************************************
*          general time loop            *
*****************************************/
void callSysDynamics(Event *event,Event *mevent);
/****************************************
*          1 host time loop             *
*****************************************/
void callSysDynamics1H(Event *event);
/****************************************
*  Free allocated memory for vectors    *
*  that are part of an Event struct     *
*****************************************/
void freeVecsEvent(Event *event);
#endif
