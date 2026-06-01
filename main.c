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
Event event;//struct for birth-death events: rates, cumulative probabilities and number of events
Event mevent;//struct for host migration events
/****************Program's Routines*****************************/
int main(void){
        
	/*setting the system*/
	callSetSystem();
	/**********/
        
#ifdef MULTIPLE_COSTS_WxT
	averInvXtMultipleCosts(&event,&mevent);
#endif
#ifdef TMEAS
	allocateMemTM();
	stime->saveT=0.;
	stime->Tf=TF;
	#ifdef AVERINVxT
	int i;
	for(i=0; i<SAMPLE; ++i){
        	openFiles();
		callSysDynamics(&event,&mevent);
		closeFiles();
		stime->Tnow=0.;
		stime->saveT=0.;
		setCI();
	}
	#endif
	#ifdef SAVE_CONFIG	
	callSysDynamics(&event,&mevent);
	#endif
	#ifdef INV_DIST
	callSysDynamics(&event,&mevent);
	#endif
	#ifdef CORRxT
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif
	#ifdef CLUSTERSxT
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif
	#ifdef CLUSTERS_DISTxT
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif
	#ifdef BESTCLUSTER_TIMES
	int i;
	for(i=0; i<SAMPLE; ++i){
        	openFiles();
		callSysDynamics(&event,&mevent);
		closeFiles();
		stime->Tnow=0.;
		setCI();
	}
	#endif
	#ifdef NUMHEVENTSxT
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif
	#ifdef DIFBACOMPxT
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif
	freeMemTM();
#endif
#ifdef STEADY_STATE_MEAS
	int i;
	stime->saveT=0.;
	for(i=0; i<SAMPLE; ++i){
		#ifdef AVINVxRH
		averInvXrh(&event,&mevent);
		#endif
		#ifdef AVINVxGH
		averInvXgh(&event,&mevent);
		#endif
		#ifdef AVINVxMB
		averInvXmb(&event,&mevent,1e-6,0.004);
		#endif
		#ifdef AVINVxCOST
		averInvXcost(&event,&mevent);
		#endif
		#ifdef AVINVxMH
		averInvXmh(&event,&mevent);
		#endif
		#ifdef RHxMHxAVINV
		rhXmhXw(&event,&mevent);
		#endif
		#ifdef COSTxMBxAVINV
		costXmbXw(&event,&mevent);
		#endif
		#ifdef COSTxPLRxAVINV
		costXplrXw(&event,&mevent);
		#endif
	}
#endif

	freeMemory();
	return 0;
}
/******************************************
 *          set system                    *
 ******************************************/
void callSetSystem(void){

        setSystem();

        return;
}
/******************************************
* 	Clean Allocated Memory            *
*******************************************/
void freeMemory(void){

	if(host){
		free(host);
	}
	if(bac){
		free(bac);
	}
	#if (Tneg>0)
	if(costvec){
		free(costvec);
	}
	#endif
        
	if(dtVec){
		free(dtVec);
	}

	//Structs and their arrays
	if(spar->inv){
		free(spar->inv);
	}
	if(spar->micr){
		free(spar->micr);
	}
	if(spar){
		free(spar);
	}
	if(stime){
		free(stime);
	}

#if (NETWORK!=0)
	int i;
	for(i=0; i<N; ++i){
		if(neighbor[i]){
			free(neighbor[i]);
		}
	}

	if(neighbor){
		free(neighbor);
	}

	if(rho_e){
		free(rho_e);
	}

#endif
	if(listh->vec){
		free(listh->vec);
	}
	if(inverselisth){
		free(inverselisth);
	}

	/*structs*/
	if(listh){
		free(listh);
	}
	if(sysmeas){
		free(sysmeas);
	}

	return;
}
