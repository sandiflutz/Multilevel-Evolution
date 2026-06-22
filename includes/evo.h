#ifndef EVO_H
#define EVO_H
/****************************************************************
*     calculates the accumulated investment in host @index      *
*****************************************************************/
double calcAcumInvest(int index);
/********************************************************
*       Set host birth-death rates without any spatial  *
*       characteristics                                 *
*********************************************************/
void setPureHostBDRates(Event *event);
/********************************************************
*       Set the spatial of the host rates:              *
*       -include migration rates                        *
*       -multiply birth rates by a step funcion of the  *
*       density of empty sites in the neighborhood of   *
*       each host                                       *
*********************************************************/
void setHostRatesSpatialEv(Event *event,Event *mevent);
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
*       2 states in the network         *
*       exchange places                 *
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
/****************************************************************************************
*       Choose a site to migrate within a specific connection distance                  *
*****************************************************************************************/
int chooseMigSiteConDis(int id, int rf);
/****************************************************************************************
*               store site network id's for host migration and                          *
*               the related probability of being chosen                                 *
*               depending on dilution.                                                  *
*               Sites included are only the one at a given                              *
*               distance                                                                *
*****************************************************************************************/
double findMigSiteswithEuclRmig(int id0,int r,int i0,int imax,double da,int *which_host,double *prob);
/****************************************************************
*       Choose site for migration within a distance r:          *
*       states of sites @idm and one chosen are exchanged.      *
*       The probability of being chosen given that a site       *
*       is within a distance r depends of on vacancy:           *
*               -prob[i]=1/sumprob, if neighbor is empty,       *
*               -prob[i]=rho_e[i]/sumprob, if site i not        *
*                       empty and rho_e[i]>0 for at least       *
*                       one of the sites being consider         *
*               -prob[i]=1/num_sites, for all sites being       *
*               consider otherwise                              *
*****************************************************************/
int chooseMigSite(int idm,int rmig);
/****************************************************************
*       Long-range random host migration:                       *
*       randomly choose a site for migration under the          *
*       the following constraints:                              *
*               -its distance from the focus site has           *
*               to be at equal or larger than a minimum         *
*               value                                           *
*               -the fraction of empty sites in it              *
*               neighborhood has to be larger than 0.           *
*****************************************************************/
int longRangeMigSite(int idm,int rmig,int rmin);
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
