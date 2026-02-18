#ifndef GLOBALS_H
#define GLOBALS_H
#include<assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include"tools.h"
/************************************************************************************
*                 Defining global constants and macros                              *
*************************************************************************************/

/**paramenters to define the structure of the system*********/
#define L               100        /*linear number of sites (square lattice case*/
#define N               (L*L)      /*number of sites*/
#define TYPES           100        /*number of types of microbe*/
#define Tplus           (TYPES-1)  /*number of types of microbe*/
#define IDBH            0          /*type of the bacteria that helps other bacteria*/  
#define NETWORK         0          /*0: well-mixed
                                    *1: square-lattice*/
#if (NETWORK==0)
        #define VIZ    N	/*number of neighbors in the well-mixed*/
#elif
        #define VIZ    4	/*number of neighbors in the square-lattice: 4 or 8*/
#endif
#define CI             1	/*0: system starts with types being randomly distributed with a uniform distribution
                       		*1: system starts with types being randomly distributed using a normal distribution for the frequencies of each type
                       		*2: system starts with only the 1 host*/

/**dynamics parameters*************/
#define Mu        0.01          /*mutation rate*/
#define Theta     1e-05         /*migration rate*/
#define Beta      1.            /*birth rate for neutral bacteria*/
#define Delta     1.            /*death rate for microbes*/
#define K_bac     (Beta/Delta)  /*bacteria carrying capacity*/
#define K_H       500
#define Gamma     0.01          /*cost for helping when the investment is 1*/
#define Gh        100           /*number of microbial generations per host generation*/
#define Sb        1.            /*strength of the dependence of hosts births on their microbial content*/
#define Sd        0.            /*strength of the dependence of hosts deaths on their microbial content*/
#define Dt_ref    0.05          /*time step for bacteria evolution*/
#define Bacv      1e-03         /*initial density of vertically transmitted microbes in a new host*/
#define Bac0      1.            /*initial bacteria density in each host (t=0)*/
#define H0        K_H           /*initial number of hosts*/
#define BSAMPLES  10            /*number of bacteria samples passed from a parent host to its offspring*/
#define SIGMA     0.05          /*variance of the trucated normal distribution for the inheritance of helpful microbes*/
#define MEANinv0  0.
#define STDinv0   0.01
#define DTVSIZE   29
/****parameters for measures/sampling and related things****************/
#define TF        1000.           /*host maximum time (measured using continuous values for the times steps)*/
#define NTS       10e7            /*maximum number of timesteps*/
#define FIG_EXT   0               /*Extension of the figure files (for the routines that create gnuplot scripts):
				  * 0:png (good for creating animations later)
                                  * 1:eps*/
#define NF       100             /*number of files creater for the routines that create one file per timestep (interval between file=(Tf_me-T0_me)/NF)*/
#define NTf_me   NTS             /*time to stop a measure*/
#define NT0_me   0               /*time to start a measure*/
/***Routine Choices***********************************************/
#if (TYPES==2)
  #define INV     1         /*(don't change it)! Investment function when there is just 2 types is given by inv[j]=(j-Tmin)/Tplus):  
			     * 0 or 1: if the types are neutral and positive
			     * 0 and -1: if the types are neutral and negative*/
#else
  #define INV     0         /*investment function:
			     *0 (paper version): inv[j]=(2(j+1)-1)/(2*TYPES)
			     *1: inv[j]=(j-Tmin)/Tplus (where Tplus=# of positive types, Tmin=# of negative types=TYPES-1-Tplus)*/
#endif
#define TV        0         /*rule for vertical transmission: 0=normal dist. (around parent bac. type freq.),1=poisson dist. for the sample size */
#if  defined(DENSb1xT)||defined(AVERINVxT)||defined(SAVE_CONFIG)||defined(INV_DIST)
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
        double mig;//migration rate
	double sb;
	double sd;
	double sigma;
	double *s;
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
        double Tnow;
	double saveT;
	int NTi;
	int NTf;
	int NTnow;
        int nfiles;
        int idh_h1;
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
extern double *dtVec;
extern char *fdatapath;
extern DynList *listh;
extern Dyn2DList *listh2d;
extern DynList *alive_viz;
extern SysParams *spar;
#endif
