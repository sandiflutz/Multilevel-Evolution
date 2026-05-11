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
	stime->Tf=TF;
	#if SAVE_CONFIG	
	callSysDynamics(&event,&mevent);
	#endif
	#if INV_DIST
	callSysDynamics(&event,&mevent);
	#endif
	#if CORRxT
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif
	#if CLUSTERSxT
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif
	#if RVNxTxCORRBAC
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif
	#if NUMHEVENTSxT
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif
	#if GENTIME
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif
	#if DIFBACOMPxT
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
	#endif

	#if defined(AVERINVxT)||(AVERINVRATExT)
	int i;
	stime->Tf=TF;
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
