/*
 * This file is part of the project Multilevel Evolution.
 *
 * This program uses the GNU Scientific Library (GSL), which is licensed
 * under the GNU GPL v3 or later. Therefore, this file and any binary
 * linked with GSL must be distributed under the terms of the GPL v3+.
 *
 * Copyright (C) 2025 Alessandra F. Lütz
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


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
#define L               100       /*linear number of sites (square lattice case*/
#define N               (L*L)    /*number of sites*/
#define TYPES           100      /*number of types of microbe*/
#define Tplus           TYPES-1  /*number of types of microbe*/
#define IDBH            0        /*type of the bacteria that helps other bacteria*/  
#define NETWORK         0        /*0: well-mixed
                                  *1: square-lattice
                                  */
#if (NETWORK==0)
        #define VIZ    (N-1)     /*number of neighbors in the well-mixed*/
#elif
        #define VIZ    4         /*number of neighbors in the square-lattice: 4 or 8*/
#endif
#define CI             1         /*0: system starts with types being randomly distributed with a uniform distribution
                                  *1: system starts with types being randomly distributed using a normal distribution for the frequencies of each type
                                  *2: system starts with only the 1 host
                                  */

/**dynamics parameters*************/
#define Mu        0.01  /*mutation rate*/
#define Theta     0.00001       /*migration rate*/
#define Beta      1.            /*birth rate for neutral bacteria*/
#define Delta     1.            /*death rate for microbes*/
#define K_bac     1.            /*bacteria carrying capacity*/
#define K_H       500
#define Gamma     0.01          /*cost for helping when the investment is 1*/
#define Gh        100           /*number of microbial generations per host generation*/
#define Sb        1.            /*strength of the dependence of hosts births on their microbial content*/
#define Sd        0.            /*strength of the dependence of hosts deaths on their microbial content*/
#define Dt_ref    0.05          /*time step for bacteria evolution*/
#define Bacv      0.001         /*initial density of vertically transmitted microbes in a new host*/
#define Bac0      1.             /*initial bacteria density in each host (t=0)*/
#define BSAMPLES  10            /*number of bacteria samples passed from a parent host to its offspring*/
#define SIGMA     0.05          /*variance of the trucated normal distribution for the inheritance of helpful microbes*/
#define MEANinv0  0.
#define STDinv0   0.01
/****parameters for measures/sampling and related things****************/
#define TF        1000.          /*host maximum time (measured using continuous values for the times steps)*/
#define NTS      10000000       /*maximum number of timesteps*/
#define SAVE_CONFIG_ID     0   /*0:png
                                *1:eps
                                */
#define INT     100            /*interval between measures*/
#define Tf_me   10000          /*time to stop a measure*/
#define T0_me   0              /*time to start a measure*/
/***Routine Choices***********************************************/
#if (TYPES==2)
  #define INV     1         /*(don't change it)! Investment function when there is just 2 types is 0 or 1 (option 1 is inv[j]=j/(TYPES-1) )*/
  #else
  #define INV     0         /*investment function:
                                 *0 (paper version): inv[j]=(2(j+1)-1)/(2*TYPES)
                                 *1: inv[j]=j/(TYPES-1)
                                 */
#endif
#if  defined(DENSb1xT)||defined(AVERINVxT)||defined(SAVE_CONFIG)||defined(TIME_VARS)
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
        double Ti;
        double Tf;
        double Tnow;
	double saveT;
	int NTi;
	int NTf;
	int NTnow;
        int interv;
        int idh_h1;
	char *ftname_pars;
	int ftnpars_size;
} TimeMeasures;
/***************************************************
 *            Global Variables                     *
 ***************************************************/
extern int *s;
extern int *host;
extern int **neighbor;
extern double **bac;
extern double *micr;
extern double *inv;
extern FILE *fvarsXt;
extern FILE *fdensb1Xt;
extern FILE *finvCumul;
extern DynList *listh;
extern DynList *listnb;
extern DynList *alive_viz;
#endif
