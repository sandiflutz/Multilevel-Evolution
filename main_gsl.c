/* main_gsl.c */

#include<stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include"globals.h"
#include"randgen_ufrgs.h"
#include"init.h"
#include"bac_gsl.h"
#include"evo.h"
#include"measures.h"
/*****************declaring routines order**************************/
void freeMemory(void);
/******global variables*********************************************/
int EXIT_N=0;
Event event;
TimeMeasures meas;
gsl_odeiv2_driver *driver=NULL;
gsl_odeiv2_system sys;
/****************Program's Routines*****************************/
int main(void){
	int i,idh,numsteps,nh,nb,dnt_h,dnt_b,dnt;
	double dt,sumprobs;
        
	/*setting the system*/
	setSystem(&event,&meas);
	
        /*Open Files for time measures*/
	#ifdef TMEAS
	strncat(meas.ftname_pars, "_gil.dat",meas.ftnpars_size-strlen(meas.ftname_pars)-1);
        openFiles(&meas);
        #endif

        /*setting microbial solver*/
        bac_make_system(&sys, spar);
        driver =
                gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_msbdf,
                1e-6,// initial step size guess
                1e-6,// absolute tolerance
                1e-6// relative tolerance
        );
	/**********/
	
	numsteps=0;
	event.timeE=0.;
	while(numsteps<=meas.NTf){
		meas.NTnow=numsteps;
		meas.Tnow=event.timeE;
		nh=listh->usize;
		event.sizeE=nh*2;
		#ifdef TMEAS
		measures(meas);
		#endif
		calcHostEvents(event.ratesE);
		cumulProb(event.sizeE,event.ratesE,event.cprobE);
		
                dt=adjustTimeStep(event.cprobE[event.sizeE-1]);//adjust dt to calculate host event probabilities
                sumprobs=dt*event.cprobE[event.sizeE-1];
                event.dtE=gillespieTime(sumprobs);//use sum of event probs. to calculate gillespie time

		calcNumSteps(&dnt_h,&dnt_b,&dnt,event.dtE,Dt_ref);
		for(i=0; i<dnt_h; ++i){
			do{
				event.whichE=selectEventCP(event.cprobE,event.sizeE);
				idh=listh->vec[event.whichE%nh];
			}while(host[idh]==2);
		
			dynamicsHost(&event);
		}
		EXIT_N=bacDynamics(driver,&event);
			
		/*no newborns anymore*/
		nb=listnb->usize;
		for(i=0; i<nb; ++i){
			host[listnb->vec[i]]=1;//no newborns anymore
		}
		listnb->usize=0;

		event.timeE+=event.dtE;
		
                numsteps+=dnt;
        }


	freeMemory();
	return EXIT_N;
}
/******************************************
 * Clean Allocated Memory                 *
 *****************************************/
void freeMemory(void){
	int i;

	free(host);
	for(i=0; i<VIZ; ++i){
		free(neighbor[i]);
	}
	free(neighbor);
	for(i=0; i<N; ++i){
		free(bac[i]);
	}
	free(bac);
	free(dtVec);
	
	//Structs and their arrays
	free(spar->s);
	free(spar->inv);
	free(spar->micr);
	free(spar);
#ifdef TMEAS
	free(meas.ftname_pars);
#endif
	free(dtVec);

#if (NETWORK!=0)
	free(alive_viz->vec);
	free(alive_viz);
#endif

	free(event.ratesE);
	free(event.cprobE);
	free(listnb->vec);
	free(listh->vec);

	/*structs*/
	free(listh);
	free(listnb);
	free(driver);

	/*close files*/
#ifdef TIME_VARS
	fclose(fvarsXt);
#endif
#ifdef DENSb1xT 
	fclose(fdensb1Xt);
#endif
#ifdef AVERINVxT
	fclose(finvCumul);
#endif
	
	return;
}
