#ifndef GLOBALS_H
#define GLOBALS_H
#include <stdbool.h>
#include<assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include"tools.h"
/************************************************************************************
*                 Defining global constants and macros                              *
*************************************************************************************/
//parameters for the dynamics 
#define TAB             	2		/*for choosing which set of values for some of the parameters to use: 1 or 2 for the specified fixed values below, and anything else for manually choosing*/

#if (TAB==1)
	#define Mu        	1e-09		/*mutation rate (tab1:1e-09, tab2:1e-02)*/
	#define Theta     	1e-06		/*migration rate (tab1:1e-06, tab2:1e-05)*/
	#define K_H      	5000		/*carrying capacity for the host layer (500 for most cases, but 5000 for fig2, types=2: if K_H=5000, use L~224 to have (L^2/K_H >=10)*/
	#define Gh        	10		/*# of microbial generations per host generation (usually 100, but fig2 uses 10, for types=2)*/
	#define Bacv      	1e-04		/*density of vertically transmitted microbes in a new host (tab1:1e-04, tab2:1e-03)*/
#elif (TAB==2)
	#define Mu        	1e-02		/*mutation rate (tab1:1e-09, tab2:1e-02)*/
	#define Theta     	1e-05		/*migration rate (tab1:1e-06, tab2:1e-05)*/
	#define K_H      	500		/*carrying capacity for the host layer (500 for most cases, but 5000 for fig2, types=2: if K_H=5000, use L~224 to have (L^2/K_H >=10)*/
	#define Gh        	100		/*# of microbial generations per host generation (usually 100, but fig2 uses 10, for types=2)*/
	#define Bacv      	1e-03		/*density of vertically transmitted microbes in a new host (tab1:1e-04, tab2:1e-03)*/
#else
	#define Mu        	1e-02		/*mutation rate (tab1:1e-09, tab2:1e-02)*/
	#define Theta     	1e-05		/*migration rate (tab1:1e-06, tab2:1e-05)*/
	#define K_H      	500		/*carrying capacity for the host layer (500 for most cases, but 5000 for fig2, types=2: if K_H=5000, use L~224 to have (L^2/K_H >=10)*/
	#define Gh        	100		/*# of microbial generations per host generation (usually 100, but fig2 uses 10, for types=2)*/
	#define Bacv      	1e-03		/*density of vertically transmitted microbes in a new host (tab1:1e-04, tab2:1e-03)*/
#endif

#define SIGMA     	0.01		/*variance of the trucated normal distribution for the inheritance of helpful microbes*/
#define DTVSIZE   	19		/*number of possible time steps (<Dt_ref=microbial time step=0.05) that can be chosen for
					 *the host dynamics so the probability of 2 consecutive host events during a @Dt_ref time interval is <0.01
					 *paper uses 19 when TYPES=2 and 29 otherwise (??)*/

#define Mh        	0.		/*host migration coeficient (>=0.): if Mh=0., there is no host migration*/

/**cost**/
#define Gamma     	1e-02		/*cost for helping when the investment is 1*/
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


/***************System**Structure*********************************************************************************************************************************************************************/
/***host**layer******/
#if (TAB==1)
	#define L               224	/*linear size of the system (for the square lattice case)*/
#elif (TAB==2)
	#define L               100
#else
	#define L               100
#endif

#define N               (L*L)		/*size number of sites*/
#define NETWORK         1		/*0: well-mixed
					 *1: square-lattice*/
#if (NETWORK==0)
        #define VIZ	N		/*number of neighbors in the well-mixed*/
#else
        #define VIZ	4		/*number of neighbors in the square-lattice: 4 or 8*/
#endif
/***bacteria****/
#define TYPES           100		/*number of types of microbe*/
#define Tpos            (TYPES-1)	/*number of positive types of microbes: positively affect host reproduction success*/
#define Tneg            0		/*number of negative types of microbes: negatively affect host reproduction success*/
//initial conditions
#if ((TYPES==2)&&(Tpos==TYPES-1))
	#define CI	3
#else
	#define CI	1		/*0: system starts with types being randomly distributed with a uniform distribution
					 *1 (default): system starts with types being randomly distributed using a normal distribution for the frequencies of each type
					 *2: system starts with only the 1 host
					 *3: all types of bacteria start with a fixed fraction of 1/TYPES */
#endif
/***Routine Choices*********************************************************************************************************************************************************************************/
#if ((TYPES==2)&&(Tpos==TYPES-1))
	#define INV     	1       /*(don't change it)! Investment function when there is just 2 types is given by inv[j]=(j-Tmin)/Tpos):  
					* 0 or 1: if the types are neutral and positive
					* 0 and -1: if the types are neutral and negative*/
	#define MUT_BIRTH_DYN  	0	/*0 (my version): mutation from type 0 to type 1, or from @TYPES-1 to @TYPES-2 happen with rate mu 
					* 1 (paper version for TYPES=100): mutation from type 0 to type 1, or from @TYPES-1 to @TYPES-2 happen with rate mu/2 
					* (all other types mutate from j to j+1 with rate mu/2 and from j to j-1 with rate mu/2) 
					*/
#else
  	#define INV     	0       /*investment function:
			     		*0 (paper version): inv[j]=(2(j+1)-1)/(2*TYPES)
			     		*1: inv[j]=(j-Tneg)/Tpos (where Tpos=# of positive types, Tneg=# of negative types=TYPES-1-Tpos)*/
	#define MUT_BIRTH_DYN  	0	/*0 (my version): mutation from type 0 to type 1, or from @TYPES-1 to @TYPES-2 happen with rate mu 
					* 1 (paper version for TYPES=100): mutation from type 0 to type 1, or from @TYPES-1 to @TYPES-2 happen with rate mu/2 
					* (all other types mutate from j to j+1 with rate mu/2 and from j to j-1 with rate mu/2) 
					*/
#endif
#define TV			0	/*rule for vertical transmission: 0=normal dist. (around parent bac. type freq.),1=poisson dist. for the sample size */
/****parameters for measures/sampling and related things*************************************************************/
#define TF		250000.         /*host maximum time (measured using continuous values for the times steps)*/
#define NTS		10e7            /*maximum number of timesteps*/
#define FIG_EXT		0               /*Extension of the image files that are gonna be used in gnuplot scripts:
					 * 0:png (good for creating animations later)
					 * * 1:eps*/
#define NF		1000             /*number of files for routines that create scripts for images*/
#define NInterv		100             /*Ninterv*Dt_ref=time interval between snapshots taken*/
#define NTf_me		10000            /*time (in #of time steps) to stop a measure*/
#define NT0_me		0               /*time to start a measure*/
#define SAMPLE		100             /*number of files with raw data that are going to be produce for measurements that require it*/
/********************************************************************************************************************************************************************/
/*****Fixed Parameters**********/
#define Beta      1.            /*birth rate for neutral bacteria*/
#define Delta     1.            /*death rate for microbes*/
#define K_bac     (Beta/Delta)  /*bacteria carrying capacity*/
#define Sb        1.            /*strength of the dependence of hosts births on their microbial content*/
#define Sd        0.            /*strength of the dependence of hosts deaths on their microbial content*/
#define Dt_ref    0.05          /*time step for bacteria evolution*/
#define Bac0      1.            /*initial bacteria density in each host (t=0)*/
#define H0        K_H           /*initial number of hosts*/
#define BSAMPLES  10            /*number of bacteria samples passed from a parent host to its offspring*/
#define MEANinv0  0.
#define STDinv0   0.01
/***LABELING (nothing to change here)***/
//square lattice labels
#define UP		0 /*label of the top neighbor (for the square lattice)*/
#define RIGHT		1 /*label of the right neighbor (for the square lattice)*/
#define DOWN		2 /*label of the bottom neighbor (for the square lattice)*/
#define LEFT		3 /*label of the left neighbor (for the square lattice)*/
#define LeftUp		4 /*label of the neighbor at the left-up diagonal (for the square lattice)*/
#define RightUp		5 /*label of the neighbor at the right-up diagonal (for the square lattice)*/
#define RightDown	6 /*label of the neighbor at the right-down diagonal (for the square lattice)*/
#define LeftDown	7 /*label of the neighbor at the left-down diagonal (for the square lattice)*/
/******defining a main MACRO for measures made during time evolution******/
#if  defined(DENSb1xT)||defined(AVERINVxT)||defined(SAVE_CONFIG)||defined(INV_DIST)||defined(MEANBFRACxT)||defined(NUMHEVENTSxT)||defined(CORRxT)||defined(GENTIME)
	#define TMEAS
#endif
/********************************************
*  Struct for System Parameters             *
*********************************************/
typedef struct{
	int kh;
	int gh;
	int idhost;
	int idbac_help;
	double kbac;
        double mu;//bacteria mutation rate
        double cost;//cost of helping for an ideal helper bacteria
        double beta;//birth rate of a neutral bacteria
        double delta;//bacteria death rate
        double mig;//bacteria migration rate
        double migh;//host migration coeficient
	double sb;
	double sd;
	double sigma;
	double *inv;
	double *micr;
} SysParams;

typedef struct{
        int sizeE;
        int usizeE;
        int whichE;
        double dtE;
        double timeE;
        double *ratesE;
        double *cprobE;
} Event;
/*******************************************
*  Structs related to Measures             *
*******************************************/
typedef struct{
        double Ti;//time of the first time measure
        double Tf;//final time
        double Tnow;//current time
	double saveT;//save data at this time
	double dth;//time interval for the host layer
	double timegh;
	int ngh;
	int NTf;//final time in number of steps
	int NTnow;//current time in number of steps
	int nfiles;//number of data files already created
        int idh_h1;
	int numb;
	int numd;
	char *ftname_pars;
	int ftnpars_size;
	FILE *file_tmeas;
} TimeMeasures;
/***************************************************
 *            Global Variables                     *
 ***************************************************/
extern int *host;
extern int **neighbor;
extern int **clneighbor;
extern int *inverselisth;
extern double **bac;
extern double **costvec;
extern double *dtVec;
extern double *timeb;
extern char *fdatapath;
extern DynList *listh;
extern DynList *alive_viz;
extern SysParams *spar;
#endif
