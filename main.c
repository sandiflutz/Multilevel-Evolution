/* main.c */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"globals.h"
#include"randgen_ufrgs.h"
#include"init.h"
#include"bac.h"
#include"evo.h"
#include"measures.h"
/*****************declaring routines order**************************/
void callSetSystem(void);
void callSysDynamics(double tf);
void callSysDynamics1H(double tf);
void freeMemory(void);
/******global variables*********************************************/
Event event;
TimeMeasures meas;
/****************Program's Routines*****************************/
int main(void){
        
	/*setting the system*/
	callSetSystem();
	/**********/
        
#ifdef TMEAS
	allocateMemTM(&meas);
	meas.saveT=meas.Ti;
	meas.Ti=0.;
	#if defined(SAVE_CONFIG)||defined(INV_DIST)
	meas.Tf=NF*NInterv*Dt_ref+1.;	
	meas.nfiles=0;
	callSysDynamics(meas.Tf);
	#elif defined(CORRxT)||defined(NUMHEVENTSxT)||defined(GENTIME)||(DIFBACOMPxT)
        openFiles(&meas);
	meas.Tf=50000.;
	callSysDynamics(meas.Tf);
	closeFiles(&meas);
	#else
	int i;
	for(i=0; i<SAMPLE; ++i){
        	openFiles(&meas);
		
		#if (CI!=2)
		callSysDynamics(meas.Tf);
		#else
		callSysDynamics1H(meas.Tf);
		#endif
		
		setCI(&meas);
		closeFiles(&meas);
	}
	#endif
	freeMemTM(&meas);
#endif
#ifdef TESTE_DYN
	meas.Tf=Dt_ref;
	callSysDynamics(meas.Tf);
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
/****************************************
*          general time loop		*
*****************************************/
void callSysDynamics(double tf){
	int i,numsteps,nh,nevents,dnumsteps;

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
			#endif
		#endif
		
		calcHostEvents(&event);
		#if (EVO==0)//complete graph with adjustable host dt
		dnumsteps=hostNTSPerBacNTS(&event);
		nevents=event.usizeE;
		for(i=0; i<nevents; ++i){
			event.cprobE[i]*=event.dtE;
		}
		evolveHostCG(dnumsteps,&event,&meas);
		#elif (EVO==1)//square lattice with adjustable host dt
		dnumsteps=hostNTSPerBacNTS(&event);
		nevents=event.usizeE;
		for(i=0; i<nevents; ++i){
			event.cprobE[i]*=event.dtE;
		}
		evolveHostSL(dnumsteps,&event,&meas);
		#endif
		nh=listh->usize;

		evoBac(Dt_ref,event.timeE);
		
		event.timeE+=Dt_ref;
                numsteps+=dnumsteps;
		#ifdef TMEAS
		meas.dth=event.dtE;
		#endif
        }

        return;
}
/****************************************
*          1 host time loop 		*
*****************************************/
void callSysDynamics1H(double tf){
	int numsteps;

	numsteps=0;
	event.timeE=0.;
	while(event.timeE<=tf){
		#ifdef TMEAS
		meas.NTnow=numsteps;
		meas.Tnow=event.timeE;
		measures(&meas);
		#endif
		evoBac(Dt_ref,event.timeE);
		event.timeE+=Dt_ref;
                ++numsteps;
        }

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
        for(i=0; i<N; ++i){
		free(costvec[i]);
        }
	free(costvec);
        
	free(dtVec);
	free(fdatapath);

	//Structs and their arrays
	free(spar->inv);
	free(spar->micr);
	free(spar);

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
