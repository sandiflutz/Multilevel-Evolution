/* main.c */
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include"globals.h"
#include"randgen_ufrgs.h"
#include"init.h"
#include"bac.h"
#include"evo.h"
#include"measures.h"
/*****************declaring routines order**************************/
void callSetSystem(void);
void callSysDynamics(double tf);
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
	#if defined(SAVE_CONFIG)||defined(INV_DIST)
	meas.Tf=NF*Dt_ref+1.;	
	callSysDynamics(meas.Tf);
	#else
	for(i=0; i<SAMPLE; ++i){
		setCI(&meas);
        	openFiles(&meas);
		
		callSysDynamics(meas.Tf);
		
		closeFiles(&meas);
	}
	#endif
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
void callSysDynamics(double tf){
	int numsteps,nh,dnumsteps;
	int *inverselisth_tmp=(int *)calloc(N,sizeof(int));
	DynList listh_tmp;
	listh_tmp.vec=(int *)calloc(N,sizeof(int));
	listh_tmp.size=listh->size;

	#ifdef NUMHEVENTSxT
	meas.numb=0;
	meas.numd=0;
		#if (NETWORK==1)
			for(int i=0; i<L; ++i){
				meas.numb_lin[i]=0;
				meas.numd_lin[i]=0;
				meas.numb_col[i]=0;
				meas.numd_col[i]=0;
			}
                #endif
	#endif	

	numsteps=0;
	event.timeE=0.;
	nh=listh->usize;
	dnumsteps=1;
	while((event.timeE<=tf)&&(nh>0)){
		#ifdef TMEAS
		meas.NTnow=numsteps;
		meas.Tnow=event.timeE;
		measures(&meas);
			#ifdef NUMHEVENTSxT
			meas.numb=0;
			meas.numd=0;
				#if (NETWORK==1)
				for(int i=0; i<L; ++i){
					meas.numb_lin[i]=0;
					meas.numd_lin[i]=0;
					meas.numb_col[i]=0;
					meas.numd_col[i]=0;
				}
				#endif
			#endif
		#endif
		
		#if (CI!=2)//not the single host case
		/*host dynamics (all routines in this block are from evo.c)*/
		calcHostEvents(&event);
			#if (EVO==0)//evolution using adjustment of host time step (as in the original paper)
				dnumsteps=evolveHostDtH(&event,&listh_tmp,inverselisth_tmp,&meas);
			#elif (EVO==1)
				event.dtE=Dt_ref;
				evolveHostTLP(&event,&listh_tmp,inverselisth_tmp,&meas);
			#endif
		/**************/
		#endif
		nh=listh->usize;

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
	#ifdef NUMHEVENTSxT
		#if (NETWORK==1)
                        free(meas.numb_lin);
                        free(meas.numd_lin);
                        free(meas.numb_col);
                        free(meas.numd_col);
                #endif
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
	free(inverselisth);

	/*structs*/
	free(listh);

	return;
}
