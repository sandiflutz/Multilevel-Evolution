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
        openFiles();
	callSysDynamics(&event,&mevent);
	closeFiles();
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

	if(host!=NULL){
		free(host);
		host=NULL;
	}
	if(bac!=NULL){
		free(bac);
		bac=NULL;
	}
        for(i=0; i<N; ++i){
		if(costvec[i]!=NULL){
			free(costvec[i]);
			costvec[i]=NULL;
		}
        }
	#if (Tneg>0)
	if(costvec!=NULL){
		free(costvec);
		costvec=NULL;
	}
	#endif
        
	if(dtVec!=NULL){
		free(dtVec);
		dtVec=NULL;
	}

	//Structs and their arrays
	if(spar->inv!=NULL){
		free(spar->inv);
		spar->inv=NULL;
	}
	if(spar->micr!=NULL){
		free(spar->micr);
		spar->micr=NULL;
	}
	if(spar!=NULL){
		free(spar);
		spar=NULL;
	}
	if(stime!=NULL){
		free(stime);
		stime=NULL;
	}

#if (NETWORK!=0)
	for(i=0; i<N; ++i){
		if(neighbor[i]!=NULL){
			free(neighbor[i]);
			neighbor[i]=NULL;
		}
	}

	if(neighbor!=NULL){
		free(neighbor);
		neighbor=NULL;
	}

	if(rho_e!=NULL){
		free(rho_e);
		rho_e=NULL;
	}

#endif
	if(listh->vec!=NULL){
		free(listh->vec);
		listh->vec=NULL;
	}
	if(inverselisth!=NULL){
		free(inverselisth);
		inverselisth=NULL;
	}

	/*structs*/
	if(listh!=NULL){
		free(listh);
		listh=NULL;
	}
	if(sysmeas!=NULL){
		free(sysmeas);
		sysmeas=NULL;
	}

	return;
}
