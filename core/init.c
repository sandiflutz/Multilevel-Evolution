/* init.c */

#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
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
	memset(host,0,sizeof(int)*N);
        
        bac=(double **)calloc(N,sizeof(double *));
	for(i=0; i<N; ++i){
        	bac[i]=(double *)calloc(TYPES,sizeof(double));
	}
	
	inverselisth=(int *)calloc(N,sizeof(int));

	dtVec=(double *) calloc(DTVSIZE,sizeof(double));
	//logSpacedVec(dtVec,1e-07,Dt_ref,DTVSIZE);//DTVSIZE,dtVec[0]=1e-07 and dtVec[DTVSIZE-1]=1e-02 are values used in the paper
	logSpacedVec(dtVec,-7,-2,DTVSIZE);//DTVSIZE,dtVec[0]=1e-07 and dtVec[DTVSIZE-1]=1e-02 are values used in the paper
        
	fdatapath=(char *)malloc(sizeof(char)*50);
	sprintf(fdatapath,"data_manipulation/");
	
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
	event->ratesE=(double *)calloc(event->sizeE,sizeof(double));
	memset(event->ratesE,0.,sizeof(double)*event->sizeE);
	event->cprobE=(double *)calloc(event->sizeE,sizeof(double));
	memset(event->cprobE,0.,sizeof(double)*event->sizeE);
	event->timeE=0.;
	event->dtE=Dt_ref;

	//lists structs
        listh = malloc(sizeof(DynList));
        if (!listh) { perror("malloc"); exit(1);}
        listh->vec = (int *)calloc(N,sizeof(int));
        listh->size=N;
        listh->usize=0;

	//Time measures structs
	#ifdef TMEAS
	meas->idh_h1=0.;
	meas->nfiles=0;
	meas->NTi=NT0_me;
	meas->NTf=NTf_me;
	meas->saveT=0;
	meas->ftnpars_size=250;
	meas->ftname_pars=(char *)calloc(meas->ftnpars_size,sizeof(char));
	sprintf(meas->ftname_pars,"N%d_Ty%d_Kh%d_net%d_Bv%f_cost%0.4f_sigma%0.2f_mu%f_mig%f",N,TYPES,spar->kh,NETWORK,Bacv,spar->cost,spar->sigma,spar->mu,spar->mig);
	#endif

	#if (NETWORK!=0)//not the well-mixed/complete graph case
	/*network array: e.g. neighbor[k][idh]=idh_viz (label of the k-th neighbor of host @idh is @idh_viz)*/
	neighbor=(int **)calloc(N,sizeof(int *));
	for(i=0; i<N; ++i){
                neighbor[i]=(int *)calloc(VIZ,sizeof(int));
        }
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
	int ok=0,namelen,dnl;
	unsigned long id;

	id = (unsigned long)time(NULL);

	dnl=200;
        namelen=strlen(meas->ftname_pars)+strlen(fdatapath)+dnl;

	char *name=(char *)calloc(namelen,sizeof(char));

#ifdef DENSb1xT
	/*each execution will produce a file with a different name (with a "random" id at the end of the name)*/
	while(ok==0){
		sprintf(name,"%sdensb1Xt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
		meas->file_tmeas=fopen(name,"r");
		if(meas->file_tmeas!=NULL){
			++id;
			fclose(meas->file_tmeas);
		}else{
			ok=1;
		}
	}
        sprintf(name,"%sdensb1Xt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
        meas->file_tmeas=fopen(name,"w");
	if (meas->file_tmeas==NULL) { perror("malloc"); exit(1);}

#endif
#ifdef AVERINVxT
	while(ok==0){
        	sprintf(name,"%saverInvXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
		meas->file_tmeas=fopen(name,"r");
		if(meas->file_tmeas!=NULL){
			++id;
			fclose(meas->file_tmeas);
		}else{
			ok=1;
		}

	}
        sprintf(name,"%saverInvXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
        meas->file_tmeas=fopen(name,"w");
	if (meas->file_tmeas==NULL) { perror("malloc"); exit(1);}
#endif
	free(name);
        return;
}
/*****************************************************
*   Populates Host and Microbial layer               *
*   Bacterial layer: types of bacteria               *
*   are uniformly distributed                        *
*****************************************************/
void initialStateUniD(void){
       int i,j,nh=0,ne=0;
        double norm,p_oc;

	p_oc=(double)H0/N;

        for(i=0; i<N; ++i){
		if(FRANDOM<p_oc){
                        host[i]=1;
			listh->vec[nh]=i;
			inverselisth[i]=nh;
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
			inverselisth[i]=N-1-ne;
			++ne;
			memset(bac[i],0.,sizeof(double)*TYPES);
			spar->micr[i]=0.;
		}
        }

	#if (NETWORK==0)//complete graph
	listh->usize=nh;
	#else
	#endif

        return;
}
/**********************************************************
* 	Populates Host and Microbial layer                *
* 	Bacterial layer: frequencies of each type j       *
* 	is the probability of normal                      *
* 	distribution with x=(investment[j]-mean)/stdinv   *
***********************************************************/
void initialStateNormD(void){
       int i,j,nh=0,ne=0;
        double x,*bacinit;
        double norm,p_oc;

	bacinit=(double *)calloc(TYPES,sizeof(double));
                       
       /*setting initial bacteria dist.: the same for every host*/	
	norm=0.;
	for(j=0; j<TYPES; ++j){
		x=(spar->inv[j]-MEANinv0)/STDinv0;
		bacinit[j]=normalProb(x);
		norm+=bacinit[j];//for normalization
	}
	for(j=0; j<TYPES; ++j){
		bacinit[j]*=Bac0/norm;
	}

	p_oc=(double)H0/N;

        for(i=0; i<N; ++i){
		if(FRANDOM<p_oc){//populating the network with only ~kh hosts (which is the carrying capacity in the birth_rate/death_rate sense)
                        host[i]=1;
			listh->vec[nh]=i;//live hosts labels are stored at the beginning of list of hosts
			inverselisth[i]=nh;
			++nh;
			spar->micr[i]=Bac0;
			for(j=0; j<TYPES; ++j){
				bac[i][j]=bacinit[j];
			}
		}else{
			listh->vec[N-1-ne]=i;//empty sites are stored at the end of the list of hosts
			inverselisth[i]=N-1-ne;
			++ne;
			memset(bac[i],0.,sizeof(double)*TYPES);
			spar->micr[i]=0.;
		}
        }

	listh->usize=nh;

	free(bacinit);
        return;
}
/*****************************************************
*   Populates Host layer with a single host          *
*   Bacteria yypes are uniformly distributed         *
*****************************************************/
void initialStateSingleH(TimeMeasures *meas){
       int i;
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
                bac[meas->idh_h1][i]*=Bac0/norm;
        }
        for(i=0; i<N; ++i){
                listh->vec[i]=i;
		inverselisth[i]=i;
        }
	exchange(inverselisth,0,meas->idh_h1);
	exchange(listh->vec,0,meas->idh_h1);

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

	tymi=TYPES-1-Tplus;//# of negative types
        for(i=0; i<TYPES; ++i){
                spar->inv[i]=(double)(i-tymi)/Tplus;
		
		if(i-tymi>=0){
			spar->s[i]=1.;
		}else{
			spar->s[i]=-1.;
		}
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

	/*setting hosts network*/	
	#if(NETWORK==1)
		squareLattice(neighbor,VIZ,N);
	#endif

	/*setting initial state (alive hosts and bacteria abundances)*/
	#if (CI==0)//uniform distribution
		initialStateUniD();
	#elif(CI==1)//frequencies come from normal distribution
		initialStateNormD();
	#else//single host
		initialStateSingleH(meas);
	#endif

	return;
}
