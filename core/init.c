/* init.c */

#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include"globals.h"
#include"randgen_ufrgs.h"
#include"tools.h"
#include"init.h"
/***********************************************************
*     allocate memory for global arrays and structs        *
*     initialize system parameters                         *
************************************************************/
void allocateMemory(Event *event,TimeMeasures *meas){
        int i;

        host=(int *)calloc(N,sizeof(int));
        neighbor=(int **)calloc(VIZ,sizeof(int *));
        for(i=0; i<VIZ; ++i){
                neighbor[i]=(int *)calloc(N,sizeof(int));
        }
        bac=(double **)calloc(N,sizeof(double *));
	for(i=0; i<N; ++i){
        	bac[i]=(double *)calloc(TYPES,sizeof(double));
	}
	
	dtVec=(double *) calloc(DTVSIZE,sizeof(double));
	logSpacedVec(dtVec,1e-07,1e-02,DTVSIZE);//DTVSIZE,dtVec[0] and dtVec[DTVSIZE-1] are values used in the paper
        
	/****structs****/
        //system parameters necessary for the equations of the microbial
        spar = malloc(sizeof(SysParams));
        if (!spar) { perror("malloc"); exit(1);}
	spar->kh=K_H;
	spar->gh=Gh;
	spar->kbac=K_bac;
        spar->mu=Mu;
        spar->cost=Gamma;
        spar->beta=Beta;
        spar->delta=Delta;
        spar->mig=Theta;
	spar->sb=Sb;
	spar->sd=Sd;
	spar->sigma=SIGMA;
        spar->micr=(double *)calloc(N,sizeof(double));
        spar->s=(double *)calloc(TYPES,sizeof(double));
        spar->inv=(double *)calloc(TYPES,sizeof(double));

        //host events struct
	event->sizeE=2*N;
	event->usizeE=0;
        event->ratesE=(double *)calloc(event->sizeE,sizeof(int));
        event->cprobE=(double *)calloc(event->sizeE,sizeof(int));
	event->timeE=0.;
	event->dtE=Dt_ref;

        //lists structs
        listh = malloc(sizeof(DynList));
        if (!listh) { perror("malloc"); exit(1);}
        listh->vec = (int *)calloc(N,sizeof(int));
        listh->size=N;
        listh->usize=0;

        listnb = malloc(sizeof(DynList));
        if (!listnb) { perror("malloc"); exit(1);}
        listnb->vec = (int *)calloc(N,sizeof(int));
        listnb->size=L;
        listnb->usize=0;

	//Time measures structs
	#ifdef TMEAS
	meas->idh_h1=0.;
	meas->Ti=T0_me;
        meas->Tf=Tf_me;
        meas->interv=INT;
        meas->NTi=0;
        meas->NTf=NTS;
        meas->saveT=0;
	meas->ftnpars_size=250;
	meas->ftname_pars=(char *)calloc(meas->ftnpars_size,sizeof(char));
	sprintf(meas->ftname_pars,"N%d_Ty%d_Kh%d_net%d_Bv%f_cost%0.4f_sigma%0.2f_mu%f_mig%f",N,TYPES,spar->kh,NETWORK,Bacv,spar->cost,spar->sigma,spar->mu,spar->mig);
	#endif

	#if (NETWORK!=0)
        alive_viz = malloc(sizeof(DynList));
        if (!aliveviz) { perror("malloc"); exit(1);}
        alive_viz->vec=(int *)calloc(VIZ,sizeof(int));
        alive_viz->size=VIZ;
        alive_viz->usize=0;
        #endif

        return;
}
/******************************************
 *  Open Global Files                     *
 ******************************************/
void openFiles(TimeMeasures *meas){

#ifdef TIME_VARS
        char *name=(char *)calloc(200,sizeof(char));
        sprintf(name,"varsXt_%s",meas->ftname_pars);
        fvarsXt=fopen(name,"w");
        free(name);
#endif
#ifdef DENSb1xT
        char *name=(char *)calloc(100,sizeof(char));
        sprintf(namedXt1,"densb1Xt_%s",meas->ftname_pars);
        fdensb1Xt=fopen(namedXt1,"w");
        free(name);
#endif
#ifdef AVERINVxT
        char *name=(char *)calloc(250,sizeof(char));
        sprintf(name,"averInvXt_%s",meas->ftname_pars);
        finvCumul=fopen(name,"w");
        free(name);
#endif
        return;
}
/*****************************************************
*   Populates Host and Microbial layer               *
*   Bacterial layer: types of bacteria               *
*   are uniformly distributed                        *
*****************************************************/
void initialStateUniD(void){
       int i,j,nh=0,ne=0,kh;
        double norm,p_oc;

	kh=spar->kh;

	p_oc=(double)kh/N;


        
        for(i=0; i<N; ++i){
		if(FRANDOM<p_oc){
                        host[i]=1;
			listh->vec[nh]=i;
			++nh;
			spar->micr[i]=Bac0;
                        norm=0.;
                        for(j=0; j<TYPES; ++j){
                                bac[i][j]=FRANDOM;
                                norm+=bac[i][j];
                        }
                        for(j=0; j<TYPES; ++j){
                                bac[i][j]*=Bac0/norm;
                        }
		}else{
			listh->vec[N-1-ne]=i;//empty sites are stored at the end of the list of hosts;
			++ne;
			memset(bac[i],0.,sizeof(double)*TYPES);
			spar->micr[i]=0.;
		}
        }

	listh->usize=nh;

        return;
}
/**********************************************************
* 	Populates Host and Microbial layer                *
* 	Bacterial layer: frequencies of each type j       *
* 	is the probability of normal                      *
* 	distribution with x=(investment[j]-mean)/stdinv   *
***********************************************************/
void initialStateNormD(void){
       int i,j,nh=0,ne=0,kh;
        double x;
        double norm,p_oc;

	kh=spar->kh;
	p_oc=(double)kh/N;

        for(i=0; i<N; ++i){
		if(FRANDOM<p_oc){
                        host[i]=1;
			listh->vec[nh]=i;
			++nh;
			spar->micr[i]=Bac0;
                        norm=0.;
                        for(j=0; j<TYPES; ++j){
                                x=(spar->inv[j]-MEANinv0)/STDinv0;
                                bac[i][j]=normalProb(x);
                                norm+=bac[i][j];//for normalization
                        }
                        for(j=0; j<TYPES; ++j){
                                bac[i][j]*=Bac0/norm;
                        }
		}else{
			listh->vec[N-1-ne]=i;//empty sites are stored at the end of the list of hosts;
			++ne;
			memset(bac[i],0.,sizeof(double)*TYPES);
			spar->micr[i]=0.;
		}
        }

	listh->usize=nh;


        return;
}
/*****************************************************
*   Populates Host layer with a single host          *
*   Bacteria yypes are uniformly distributed         *
*****************************************************/
void initialStateSingleH(TimeMeasures *meas){
       int i,old;
        double norm;

        
	for(i=0; i<N; ++i){
		memset(bac[i],0.,sizeof(double)*TYPES);
		spar->micr[i]=0.;
	}
	
	meas->idh_h1=0;
	host[meas->idh_h1]=1;
	spar->micr[meas->idh_h1]=Bac0;
	norm=0.;
        for(i=0; i<TYPES; ++i){
                bac[meas->idh_h1][i]=FRANDOM;
                norm+=bac[meas->idh_h1][i];
        }
        for(i=0; i<TYPES; ++i){
                bac[meas->idh_h1][i]=bac[meas->idh_h1][i]*Bac0/norm;
        }
        for(i=0; i<N; ++i){
                listh->vec[i]=i;
        }
	old=listh->vec[0];
	listh->vec[0]=meas->idh_h1;
	listh->vec[meas->idh_h1]=old;

	listh->usize=1;

        return;
}
/*************************************************
* Set Investments for a System with only helpers *
* (different types of helpers,  0<investment<1   *
* (paper version)                                *
**************************************************/
void setInvestmentsPaper(void){
        int i;

        for(i=0; i<TYPES; ++i){
                spar->inv[i]=(2.*(i+1)-1)/(2.*TYPES);
		spar->s[i]=1.;
        }

        return;
}
/*************************************************
* Set Investments (-wmin<=investment<=1)         *
**************************************************/
void setInvestments(void){
        int i,tymi;

	memset(spar->s,1.,sizeof(double)*TYPES);
	tymi=TYPES-1-Tplus;
        for(i=0; i<TYPES; ++i){
                spar->inv[i]=(double)(i-tymi)/Tplus;
		if(i-tymi<0)spar->s[i]=-1;
        }

        return;
}
/****************************************************************************
 *                     Build System                                         *
 ***************************************************************************/
void setSystem(Event *event,TimeMeasures *meas){
	
        (void)start_randomic(0);

	/*Allocating memory*/
	allocateMemory(event,meas);
	/*setting investment vector*/
	#if (INV==0) 
		setInvestmentsPaper();
	#else
		setInvestments();
	#endif
	
	#if(NETWORK==0)
		setCompleteGraph(neighbor,N);
	#elif(NETWORK==1)
		squareLattice(neighbor,VIZ,N);
	#endif

	/*initial state*/
	#if (CI==0)//uniform distribution
		initialStateUniD();
	#elif(CI==1)//frequencies come from normal distribution
		initialStateNormD();
	#else//single host
		initialStateSingleH();
	#endif

	return;
}
