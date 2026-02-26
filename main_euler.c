/* main_euler.c */

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include"globals.h"
#include"randgen_ufrgs.h"
#include"init.h"
#include"bac_euler.h"
#include"evo.h"
#include"measures.h"
/*****************declaring routines order**************************/
void callSetSystem(void);
void callSysDynamics(int nst);
void freeMemory(void);
/******global variables*********************************************/
Event event;
TimeMeasures meas;
/****************Program's Routines*****************************/
int main(void){
	int i;
        
	/*setting the system*/
	callSetSystem();
	/**********/
        
#ifdef TMEAS
	#if (EVO==0)
        strncat(meas.ftname_pars, "_eu",meas.ftnpars_size-strlen(meas.ftname_pars)-1);
	#else
        strncat(meas.ftname_pars, "_tlp",meas.ftnpars_size-strlen(meas.ftname_pars)-1);
	#endif
	for(i=0; i<SAMPLE; ++i){
        	openFiles(&meas);
		callSysDynamics(NTf_me);
		closeFiles(&meas);
	}
#endif


	freeMemory();
	return 0;
}
/******************************************
 *          set system                    *
 ******************************************/
void callSetSystem(void){

        setSystem(&event,&meas);

        return;
}
/******************************************
*          time loop                     *
******************************************/
void callSysDynamics(int nst){
	int i,idh,numsteps,nh,dnumsteps;
	double dt,sumprobs;
	int *inverselisth_tmp=(int *)calloc(N,sizeof(int));
	DynList listh_tmp;
	listh_tmp.vec=(int *)calloc(N,sizeof(int));
	listh_tmp.size=listh->size;

	#ifdef NUMHEVENTSxT
	meas.numb=0;
	meas.numd=0;
	#endif	

	numsteps=0;
	event.timeE=0.;
	nh=listh->usize;
	while((numsteps<=NTS)&&(nh>0)){
		#ifdef TMEAS
		meas.NTnow=numsteps;
		meas.Tnow=event.timeE;
		measures(&meas);
			#ifdef NUMHEVENTSxT
			meas.numb=0;
			meas.numd=0;
			#endif
		#endif
		
		#if (CI!=2)//not the single host case
		/*host dynamics (all routines in this block are from evo.c)*/
		calcHostEvents(&event);
			#if (EVO==0)//evolution using adjustment of host time step (as in the original paper)
				dnumsteps=evolveHostDtH(&event,&listh_tmp,inverselisth_tmp,&meas);
			#else
				evolveHostTLP(&event,&listh_tmp,inverselisth_tmp,&meas);
				dnumsteps=1;
			#endif
		/**************/
		#endif

		bac_euler(Dt_ref,spar);//from bac_eu.c
		
		event.timeE+=Dt_ref;
                numsteps+=dnumsteps;
        }

	free(inverselisth_tmp);
	free(listh_tmp.vec);
        return;
}
/******************************************
* 	Clean Allocated Memory            *
*******************************************/
void freeMemory(void){
	int i;

	free(host);

	for(i=0; i<N; ++i){
		free(bac[i]);
	}
	free(bac);
        
	free(dtVec);
	free(fdatapath);

	//Structs and their arrays
	free(spar->s);
	free(spar->inv);
	free(spar->micr);
	free(spar);
#ifdef TMEAS
        free(meas.ftname_pars);
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

	return;
}
