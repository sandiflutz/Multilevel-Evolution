/* main_gsl.c */

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include"globals.h"
#include"randgen_ufrgs.h"
#include"init.h"
#include"bac_gsl.h"
#include"evo.h"
#include"measures.h"
/*****************declaring routines order**************************/
void callSetSystem(void);
void callSysDynamics(int nst);
void freeMemory(void);
/******global variables*********************************************/
int EXIT_N=0;
Event event;
TimeMeasures meas;
gsl_odeiv2_driver *driver=NULL;
gsl_odeiv2_system sys;
/****************Program's Routines*****************************/
int main(void){
	int i;	
       
	/*setting the system*/
	callSetSystem();
	/**********/

#ifdef TST
	for(i=0; i<DTVSIZE; ++i){
		printf("%d %.16f\n",i,dtVec[i]);
	}
#endif
#ifdef TMEAS
	for(i=0; i<SAMPLE; ++i){
		callSysDynamics(NTf_me);
	}
#endif

	freeMemory();
	return EXIT_N;
}
/******************************************
 *  		set system                *
 ******************************************/
void callSetSystem(void){
	
	/* setting the system: setSystem() is in init.c */
	
	setSystem(&event,&meas);
	
        /*Open Files for time measures: openFiles() is in init.c*/
	#ifdef TMEAS
	strncat(meas.ftname_pars, "_gil",meas.ftnpars_size-strlen(meas.ftname_pars)-1);//add "_gil" at the end of the file name
        openFiles(&meas);
        #endif

        /*setting microbial solver (from GSL)*/
        bac_make_system(&sys, spar);
        driver =
                gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_msbdf,
                1e-6,// initial step size guess
                1e-6,// absolute tolerance
                1e-6// relative tolerance
        );
	/**********/

	return;
}
/******************************************
 *          time loop                     *
 ******************************************/
void callSysDynamics(int nst){
	int i,idh,id_list,numsteps,nh,dnt_h,dnt;
	double sumprobs;
	int *inverselisth_tmp=(int *)calloc(N,sizeof(int));
        DynList listh_tmp;
        listh_tmp.vec=(int *)calloc(N,sizeof(int));
        listh_tmp.size=listh->size;
	
	nh=listh->usize;
	
	numsteps=0;
	event.timeE=0.;
	dnt=1;
	while((numsteps<=nst)&&(nh>0)){
		#ifdef TMEAS
		meas.NTnow=numsteps;
		meas.Tnow=event.timeE;
		measures(&meas);
		#endif
		
		#if (CI!=2)//not the single host case
		/*host dynamics (all routines in this block are from evo.c)*/
		calcHostEvents(&event);/*calculates host events rates and stored them in the vector @ratesE from the struct @event
					*The corresponding cumulative probability is stored in the vector @cprobE from the same struct*/

                sumprobs=event.cprobE[event.usizeE-1];
                event.dtE=gillespieTime(sumprobs);//use sum of event probs. to calculate gillespie time

		dnt_h=ceil(Dt_ref/event.dtE);
		printf("dtE=%.7f dnt_h=%d\n",event.dtE,dnt_h);

		/*temporary host list (to keep track of the changes in the host list, that have to be updated after the host time substeps)*/
		listh_tmp.usize=listh->usize;
		for(i=0; i<N; ++i){
                        listh_tmp.vec[i]=listh->vec[i];
                        inverselisth_tmp[i]=inverselisth[i];
                }
		/******/

		for(i=0; i<dnt_h; ++i){
			event.whichE=selectEventCP(event.cprobE,event.usizeE);
			id_list=event.whichE%nh;
			idh=listh->vec[id_list];

			if(listh_tmp.usize>1){
				if(host[idh]==1){//if host is alive
					switch(event.whichE/nh){
						case 0:	dynamicsHost(0,idh,&listh_tmp,inverselisth_tmp);//reproduction	
							break;
						case 1: dynamicsHost(1,idh,&listh_tmp,inverselisth_tmp);//death
							break;
					}
				}
			}else{
				exit(1);
			}
		}
		//updating hosts dynamic list
                for(i=0; i<N; ++i){
                        listh->vec[i]=listh_tmp.vec[i];
                        inverselisth[i]=inverselisth_tmp[i];
                }
                listh->usize=listh_tmp.usize;

		nh=listh->usize;
		/*******************/
		#endif

		event.dtE=ceil(event.dtE/Dt_ref)*Dt_ref;//time ¨step¨ (how much time is going to pass, not necessarily the one used for integration) for the microbial dynamics is always >= Dt_ref 
		printf("time=%f dtE_atualizado=%.16f \n\n",event.timeE,event.dtE);
	
		bacDynamics(driver,&event);//from bac_gsl.c

		event.timeE+=event.dtE;
                numsteps+=dnt;
        }

	free(inverselisth_tmp);
	free(listh_tmp.vec);
	return;
}
/******************************************
 * Clean Allocated Memory                 *
 *****************************************/
void freeMemory(void){
	int i;

	free(host);
	for(i=0; i<N; ++i){
		free(bac[i]);
	}
	free(bac);

	free(inverselisth);

	free(dtVec);
	free(fdatapath);
	
	//Structs and their arrays
	free(spar->s);
	free(spar->inv);
	free(spar->micr);
	free(spar);
#ifdef TMEAS
	free(meas.ftname_pars);
	#if  defined(DENSb1xT)||defined(AVERINVxT)
	fclose(meas.file_tmeas);
	#endif
#endif

#if (NETWORK!=0)
	for(i=0; i<N; ++i){
		free(neighbor[i]);
	}
	free(neighbor);
	free(alive_viz->vec);
	free(alive_viz);
#endif

	free(event.ratesE);
	free(event.cprobE);
	free(listh->vec);

	/*structs*/
	free(listh);
	
	if (driver!=NULL) {
		gsl_odeiv2_driver_free(driver);
		driver=NULL;
	}
	
	return;
}
