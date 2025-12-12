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
SysParams spar;
TimeMeasures meas;
gsl_odeiv2_driver *driver=NULL;
gsl_odeiv2_system sys;
/****************Program's Routines*****************************/
int main(void){
	int i,idh,numsteps,nh,nb;
	double dt,sumprobs,time,t_tmp;
        
	/*setting the system*/
	setSystem(&spar,&event,&meas);
	
        /*Open Files for time measures*/
	#ifdef TMEAS
	strncat(meas.ftname_pars, "_gil.dat",meas.ftnpars_size-strlen(meas.ftname_pars)-1);
        openFiles(&spar,&meas);
        #endif

        /*setting microbial solver*/
        bac_make_system(&sys, &spar);
        driver =
                gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_msbdf,
                1e-6,// initial step size guess
                1e-6,// absolute tolerance
                1e-6// relative tolerance
        );
	/**********/
	
	dt=Dt_ref;
	numsteps=0;
	time=0.;
	t_tmp=0.;
	meas.Ti=0.;
	meas.NTi=0;
	meas.NTf=NTS;
	while(numsteps<=meas.NTf){
		meas.NTnow=numsteps;
		meas.Tnow=time;
		nh=listh->usize;
		event.timeE=time;
		event.sizeE=nh*2;
		#ifdef TMEAS
		measures(meas,spar);
		#endif
		calcHostEvents(event.ratesE,&spar);
		cumulProb(event.sizeE,event.ratesE,event.cprobE);
		
                dt=adjustTimeStep(event.cprobE[event.sizeE-1]);
                sumprobs=dt*event.cprobE[event.sizeE-1];
                event.dtE=gillespieTime(sumprobs);

                do{
			event.whichE=selectEventCP(event.cprobE,event.sizeE);
			idh=listh->vec[event.whichE%nh];
                }while(host[idh]==2);
		
		dynamicsHost(&event,&spar);

		EXIT_N=bacDynamics(driver,&spar,&event);

		t_tmp+=event.dtE;
		if(t_tmp>=Dt_ref){
			/*no newborns anymore*/
                	nb=listnb->usize;
               	 	for(i=0; i<nb; ++i){
                        	host[listnb->vec[i]]=1;//no newborns anymore
                	}
                	listnb->usize=0;
			t_tmp=0.;
		}

		time=event.timeE;
                ++numsteps;
        }


	freeMemory();
	return EXIT_N;
}
/******************************************
 * Clean Allocated Memory                 *
 *****************************************/
void freeMemory(void){
	int i;

	free(s);
	free(host);
	for(i=0; i<VIZ; ++i){
		free(neighbor[i]);
	}
	free(neighbor);
	for(i=0; i<N; ++i){
		free(bac[i]);
	}
	free(bac);
	free(micr);
	free(inv);
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
