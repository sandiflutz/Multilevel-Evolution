#ifndef GLOBALS_H
#define GLOBALS_H
#include<stdbool.h>
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"tools.h"
/************************************************************************************
*                 Defining global constants and macros                              *
*************************************************************************************/
#define NETWORK			1			/*0: complete graph 
							 *1: square lattice
							 */ 
#define Mh        		4.			/*host migration coeficient (>=0.): if Mh=0., there is no host migration*/

#define GR_CORR			0			/*0: effective rates are individual
							 *1: the effective host rates (for the square lattice) of a host are averaged over its focus group*/


#define L               	100			/*linear system size*/
#define N               	(L*L)			/*size number of sites*/
#define K_H      		500			/*carrying capacity for the host layer*/
#define Gh        		100			/*# of microbial generations per host generation (usually 100, but fig2 uses 10, for types=2)*/
#if (NETWORK==0)
        #define VIZ		N			/*number of neighbors in the well-mixed*/
#else
        #define VIZ		4			/*number of neighbors in the square-lattice: 4 or 8*/
#endif

/**microbial layer***/

#define TYPES           	101			/*number of types of microbes ()*/
#define Tpos            	(TYPES-1)		/*number of positive types of microbes: positively affect host reproduction success*/
#define Tneg            	0			/*number of negative types of microbes: negatively affect host reproduction success*/

#define Mu0			1e-2			/*mutation rate for TYPES=101 (if TYPES!=101, mutation rate is calculated based on this value)*/

#if (TYPES!=101)
	#define Mu        	(Mu0*(TYPES-1.0)/100.0)	/*mutation rate (default for TYPES=101 is Mu=Mu0=1e-02)*/
#else
	#define Mu		Mu0
#endif
#define Theta     		1e-06			/*migration rate (tab1:1e-06, tab2:1e-05)*/
#define Bacv      		1e-03			/*density of vertically transmitted microbes in a new host (tab1:1e-04, tab2:1e-03)*/

#define SIGMA     		0.05			/*standart deviation of the trucated normal distribution for the inheritance of helpful microbes*/
#define DTVSIZE   		19			/*number of possible time steps (<Dt_ref=microbial time step=0.05) that can be chosen for
						 	*paper uses 19 when TYPES=2 and 29 otherwise (??)*/

/*****COST*******************************************************/
#define Gamma     		0.05		/*cost for helping when the investment is 1*/
/***Parameters for cases where there are negative types*******************/
/*when there are negative types, the cost for positive and negative types can be influenced by the total frequency of the negative types (f⁻): 
 * for negative types: cost Gamma*investiment[type] is multiplied by CRnn0*exp(-CRnn1*d⁻/(1-f⁻)) 
 * for positive types: cost Gamma*investiment[type] is multiplied by CRnp0*exp(-CRnp1*f⁻/(1-f⁻)) */
//cost function for negative types
#define CRnn0           0.		/*0.: negative types have no cost 
					 *1.: cost goes to Gamma*investment[type] when f⁻->0
					 *note: investment is negative for negative types and the birth rate has this term (1-cost) which is =(1.+|cost|) for negative types, unless CRnn0<0.
					 *->3 general possibilities: 1)CRnn0>0: negatives reproduce faster than neutrals
					 * 2) CRnn0=0: negative types reproduce at the same rate as neutrals 
					 * 3) CRnn0<0: negative types reproduce slower than neutral*/ 
#define CRnn1           1.		/*for 0., the cost for negative types doesn't depend on f⁻ (it's= Gamma*investment[type]*CRnn0) */ 
//cost reduction for positive types
#define CRnp0           1.		/*0.: positive types have no cost 
					 *1.: cost goes to Gamma*investment[type] when f⁻->0*/ 
#define CRnp1           1.5		/*for 0., the cost for positive types doesn't depend on f⁻ (it's= Gamma*investment[type]*CRnp0) */ 


/***initial conditions********************************************************/
#define Bac0      		1.		/*initial bacteria density in each host (t=0)*/
#define H0        		K_H		/*initial number of hosts*/
#if ((TYPES==2)&&(Tpos==TYPES-1))
	#define CI		0
#else
	#define CI		1		/*0: system starts with types being randomly distributed with a uniform distribution
						 *1 (default): system starts with types being randomly distributed using a normal distribution for the frequencies of each type
						 *2: system starts with only the 1 host*/
#endif
/***Routine Choices*********************************************************************************************************************************************************************************/
#define TV			1	/*rule for vertical transmission:
					 *0=normal dist. (around parent bac. type freq.) 
					 *1=poisson distribuition for the number of times a type of bacteria from the parent host is chosen for the sample passed to the offspring*/
#define OFFCOMP			1	/*When DIFBACOMPxT is active: choose what to measure (related to offspring microbial composition) 
					 *0: measure of parent-offspring mean diff. in microbial composition (sample comes from the last @SAMPLE reproductions)
					 *1: measure of mean offspring accumulated investment (sample comes from the last @SAMPLE reproductions)*/
/****parameters for measures/sampling and related things*************************************************************/
#define TF			30000		/*host maximum time (measured using continuous values for the times steps)*/
#define Ttrans			15000		/*transient time (to a first trial)*/
#define	Twin			25000		/*time window for measures*/
#define FIG_EXT			0               /*Extension of the image files that are gonna be used in gnuplot scripts:
					 	* 0:png (good for creating animations later)
						* 1:eps*/
#define NF			5000		/*number of files for routines that create scripts for images*/
#define NInterv			1000		/*Ninterv*Dt_ref=time interval between snapshots taken*/
#define SAMPLE			100             /*general sample size of measures done within the program (during evolution or number of files produced with raw data)*/
#define EPS			1e-8
#define BestWtr			0.98		/*best investment threshold (used for keeping track of clusters with high investments)*/
/********************************************************************************************************************************************************************/
/*****Fixed Parameters**********/
#define Beta      		1.		/*birth rate for neutral bacteria*/
#define Delta     		1.		/*death rate for microbes*/
#define K_bac     		(Beta/Delta)	/*bacteria carrying capacity*/
#define Sb        		1.		/*strength of the dependence of hosts births on their microbial content*/
#define Sd        		0.		/*strength of the dependence of hosts deaths on their microbial content*/
#define Dt_ref    		0.05		/*time step for bacteria evolution*/
#define BSAMPLES  		10		/*number of bacteria samples passed from a parent host to its offspring*/
#define MEANinv0  		0.		/**/
#define STDinv0   		0.01		/**/
/***LABELING (nothing to change here)***/
//square lattice labels
#define UP			0 /*label of the top neighbor (for the square lattice)*/
#define RIGHT			1 /*label of the right neighbor (for the square lattice)*/
#define DOWN			2 /*label of the bottom neighbor (for the square lattice)*/
#define LEFT			3 /*label of the left neighbor (for the square lattice)*/
#define LeftUp			4 /*label of the neighbor at the left-up diagonal (for the square lattice)*/
#define RightUp			5 /*label of the neighbor at the right-up diagonal (for the square lattice)*/
#define RightDown		6 /*label of the neighbor at the right-down diagonal (for the square lattice)*/
#define LeftDown		7 /*label of the neighbor at the left-down diagonal (for the square lattice)*/
/******defining a main MACRO for measures made during time evolution******/
#if  defined(AVERINVxT)||defined(SAVE_CONFIG)||defined(INV_DIST)||defined(NUMHEVENTSxT)||defined(CORRxT)||defined(DIFBACOMPxT)||defined(CLUSTERSxT)||defined(CLUSTERS_DISTxT)||defined(BESTCLUSTER_TIMES)
	#define TMEAS
#endif
#if defined(AVINVxRH)||defined(AVINVxGH)||defined(AVINVxMB)||defined(AVINVxCOST)||defined(AVINVxMH)||defined(RHxMHxAVINV)||defined(COSTxMBxAVINV)
	#define STEADY_STATE_MEAS
#endif
/********************************************
*  Struct for System Parameters             *
*********************************************/
typedef struct{
	int kh;
	int gh;
	double kbac;
        double mu;//bacteria mutation rate
        double cost;//cost of helping for an ideal helper bacteria
        double mig;//bacteria migration rate
        double mh;//host migration coeficient
	double sigma;
	double fmin;
	double *inv;
	double *micr;
} SysParams;

typedef struct{
        int sizeE;//maximum number of host events: birth and death to each host (=2*#ofsites in the system)
        double *ratesE;//vector for the host event rates
        double *cprobE;//vector for the cumulative probabilities of each host event
} Event;
typedef struct{
	double dth;//host current timestep
	double saveT;//next time to measure something
	double tinterval;//next time to measure something
        double transtime;//transient time
	double timewindow;//time window to measure something
        double Tf;//final time
        double Tnow;//current time
}SysTimes;
/*******************************************
*  Structs related to Measures             *
*******************************************/
typedef struct{
	int numb;
	int numd;
} EvMeasures;
typedef struct{
	int nw;
	double averw;
	double averw2;
} SysMeasures;
typedef struct{
	int fnsize;
	char *fname;
	char *fdatapath;
	FILE *file;
} GenFile;
typedef struct{
	int sizeCLF;
	int whichLBF;
	double stateCL;
} ClusterFullID;//identity info of a cluster
/***************************************************
 *            Global Variables                     *
 ***************************************************/
extern int *host;
extern int **neighbor;
extern int *inverselisth;
extern int *lb;//vector of cluster labels
extern ClusterFullID *maxclw;//struct for the identity of the cluster with the best investment
extern double *rho_e;
extern double *bac;
extern double **costvec;
extern double *dtVec;
extern DynList *listh;
extern DynVec *offcomp;
extern DynVec *avinv;
extern SysParams *spar;
extern SysTimes *stime;
extern GenFile *gfile;
extern EvMeasures *meas;
extern SysMeasures *sysmeas;
extern DynList *bestwlisth;
#endif
