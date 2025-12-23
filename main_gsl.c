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
        
	/*setting the system*/
	callSetSystem();
	/**********/
	
	#ifdef TMEAS
	callSysDynamics(NTf_me);
        #endif

	freeMemory();
	return EXIT_N;
}
/******************************************
 *          set system                    *
 ******************************************/
void callSetSystem(void){
	
	/*setting the system*/
	setSystem(&event,&meas);
	
        /*Open Files for time measures*/
	#ifdef TMEAS
	strncat(meas.ftname_pars, "_gil",meas.ftnpars_size-strlen(meas.ftname_pars)-1);
        openFiles(&meas);
        #endif

        /*setting microbial solver*/
        bac_make_system(&sys, spar);
        driver =
                gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_msbdf,
                Dt_ref,// initial step size guess
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
	int i,idh,id_list,numsteps,nh,dnt_h,dnt_b,dnt;
	double dt,sumprobs;
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
		/*host dynamics*/
		calcHostEvents(&event);
		
                dt=adjustTimeStep(event.cprobE[event.usizeE-1]);//adjust dt to calculate host event probabilities
		calcNumSteps(&dnt_h,&dnt_b,&dnt,dt,Dt_ref);
		dt=Dt_ref/(double)dnt_h;
                sumprobs=dt*event.cprobE[event.usizeE-1];
                event.dtE=gillespieTime(sumprobs);//use sum of event probs. to calculate gillespie time
		calcNumSteps(&dnt_h,&dnt_b,&dnt,event.dtE,Dt_ref);
		for(i=0; i<N; ++i){
                        listh_tmp.vec[i]=listh->vec[i];
                        inverselisth_tmp[i]=inverselisth[i];
                }
		for(i=0; i<dnt_h; ++i){
			event.whichE=selectEventCP(event.cprobE,event.usizeE);
			id_list=event.whichE%nh;
			idh=listh->vec[id_list];
		
			if(host[idh]==1){
				switch(event.whichE/nh){
					case 0:	dynamicsHost(0,idh,&listh_tmp,inverselisth_tmp);	
						break;
					case 1: dynamicsHost(1,idh,&listh_tmp,inverselisth_tmp);
						break;
				}
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
	
		EXIT_N=bacDynamics(driver,&event);

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
