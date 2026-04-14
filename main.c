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
        
#ifdef TMEAS
	allocateMemTM();
	stime->saveT=0.;
	#if defined(SAVE_CONFIG)||defined(INV_DIST)	
	callSysDynamics(&event,&mevent);
	#elif defined(CORRxT)||defined(NUMHEVENTSxT)||defined(GENTIME)||(DIFBACOMPxT)
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#else
	int i;
	stime->Tf=TF;
	stime->tinterval=10.;
	for(i=0; i<SAMPLE; ++i){
		stime->saveT=stime->timewindow;
		stime->Tnow=0.;
        	openFiles();
		
		#if (CI!=2)
		callSysDynamics(&event,&mevent);
		#else
		callSysDynamics1H(&event);
		#endif
		
		fprintf(gfile->file,"\n");
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
		#ifdef AVINVxRH
		averInvXrh(&event,&mevent);
		#endif
		#ifdef AVINVxGH
		averInvXgh(&event,&mevent);
		#endif
		#ifdef AVINVxMB
		averInvXmb(&event,&mevent,1e-6,1e-2);
		#endif
		#ifdef AVINVxCOST
		averInvXcost(&event,&mevent);
		#endif
		#ifdef RHxMHxAVINV
		rhXmhXw(&event,&mevent);
		#endif
		#ifdef COSTxMBxAVINV
		costXmbXw(&event,&mevent);
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
	int i;

	free(host);

	for(i=0; i<N; ++i){
		free(bac[i]);
	}
	free(bac);
        for(i=0; i<N; ++i){
		free(costvec[i]);
        }
	#if (Tneg>0)
	free(costvec);
	#endif
        
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

	free(rho_e);

#endif
	free(listh->vec);
	free(inverselisth);

	/*structs*/
	free(listh);
	free(sysmeas);

	return;
}
