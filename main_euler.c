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
int EXIT_N=0;
Event event;
TimeMeasures meas;
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
        strncat(meas.ftname_pars, "_eu",meas.ftnpars_size-strlen(meas.ftname_pars)-1);
        openFiles(&meas);
        #endif

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

	dt=Dt_ref;
	numsteps=0;
	dnumsteps=1;
	event.timeE=0.;
	nh=listh->usize;
	while((numsteps<=NTS)&&(nh>0)){
		#ifdef TMEAS
		meas.NTnow=numsteps;
		meas.Tnow=event.timeE;
		measures(&meas);
		#endif
		
		#if (CI!=2)//not the single host case
		/*host dynamics*/
		calcHostEvents(&event);

		dt=adjustTimeStep(event.cprobE[event.usizeE-1]);
		dnumsteps=ceil(Dt_ref/dt);
		dt=Dt_ref/(double)dnumsteps;
		sumprobs=dt*event.cprobE[event.usizeE-1];
		printf("ns=%d sumprobs=%f\n",numsteps,sumprobs);
                event.dtE=dt;
		/*temporary host list (to keep track of the changes in the real time changes in the host list, that have to be update after host time substeps)*/
                listh_tmp.usize=listh->usize;
                for(i=0; i<N; ++i){
                        listh_tmp.vec[i]=listh->vec[i];
                        inverselisth_tmp[i]=inverselisth[i];
                }
                /******/
		for(i=0; i<dnumsteps; ++i){
			event.whichE=selectEventCP(event.cprobE,event.usizeE);
			idh=listh->vec[event.whichE%nh];
                	
			if((FRANDOM<sumprobs)&&(host[idh]==1)){
                                switch(event.whichE/nh){
                                        case 0: dynamicsHost(0,idh,&listh_tmp,inverselisth_tmp);
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
		/**************/
		#endif

		bac_euler(Dt_ref,spar);
		
		event.timeE+=Dt_ref;
                numsteps+=dnumsteps;
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

	return;
}
