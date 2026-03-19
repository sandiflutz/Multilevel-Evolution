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
void freeMemory(void);
/******global variables*********************************************/
Event event;
/****************Program's Routines*****************************/
int main(void){
        
	/*setting the system*/
	callSetSystem();
	/**********/
        
#ifdef TMEAS
	allocateMemTM();
	stime->saveT=0.;
	stime->Tf=stime->save+stime->timewindow+1.;	
	#if defined(SAVE_CONFIG)||defined(INV_DIST)
	stime->Tf=stime->save+stime->timewindow+1.;	
	callSysDynamics(&event);
	#elif defined(CORRxT)||defined(NUMHEVENTSxT)||defined(GENTIME)||(DIFBACOMPxT)
        openFiles();
	callSysDynamics(&event);
	closeFiles();
	#else
	int i;
	for(i=0; i<SAMPLE; ++i){
        	openFiles();
		
		#if (CI!=2)
		callSysDynamics(&event);
		#else
		callSysDynamics1H(&event);
		#endif
		
		setCI();
		closeFiles();
	}
	#endif
	freeMemTM();
#endif
#ifdef STEADY_STATE_MEAS
	int i;
	stime->saveT=0.;
	for(i=0; i<SAMPLE; ++i){
		averInvXrh(&event);
	}
#endif

	freeMemory();
	return 0;
}
/******************************************
 *          set system                    *
 ******************************************/
void callSetSystem(void){

        setSystem(&event);

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

	//Structs and their arrays
	free(spar->inv);
	free(spar->micr);
	free(spar);
	free(stime);

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
