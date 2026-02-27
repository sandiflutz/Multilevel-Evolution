/* init.c */

#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
	int e0,ef,m0,mf;//variable related to the range of dtVec: [m0*10^e0;mf*10^ef]

        host=(int *)calloc(N,sizeof(int));
	memset(host,0,sizeof(int)*N);
        
        bac=(double **)calloc(N,sizeof(double *));
	for(i=0; i<N; ++i){
        	bac[i]=(double *)calloc(TYPES,sizeof(double));
	}
	
	inverselisth=(int *)calloc(N,sizeof(int));

	dtVec=(double *) calloc(DTVSIZE,sizeof(double));
	m0=1;
	e0=-7;
	mf=5;
	ef=-2;
	logSpacedVec(dtVec,e0,ef,m0,mf,DTVSIZE);/*filling dt vector @dtVec with values from m0*10^(e0) to mf*10^(ef) logarithmic spaced
					       *passing: 1)dt vector,2)e0,3)ef,4)m0,5)mf,6)size of the dt vector  
					       *range: dtVec[0]=1e-07 and dtVec[DTVSIZE-1]=1e-02 are values used in the paper
					       *another option: dtVec[DTVSIZE-1]=0.05 is the value of the microbial time step in the original paper
					       */
        
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
	meas->Ti=0.;
	meas->Tf=TF;
	meas->NTf=NTf_me;
	meas->saveT=0;
	meas->dth=Dt_ref;
	meas->ftnpars_size=400;
	meas->ftname_pars=(char *)calloc(meas->ftnpars_size,sizeof(char));

	char *ngeral=(char *)calloc(200,sizeof(char));
	char *nevo = (char *)calloc(50,sizeof(char));
	char *ntneg = (char *)calloc(50,sizeof(char));
	char nparam[4][10];
	double param[4];
	param[0]=Bacv;
	param[1]=spar->cost;
	param[2]=spar->mu;
	param[3]=spar->mig;
	for(i=0; i<4; ++i){
		if(param[i]==0){
			sprintf(nparam[i],"0");
		}else{
			sprintf(nparam[i],"1e%d",(int)log10(param[i]));
		}	
	}
	sprintf(ngeral,"N%d_Ty%d_Tp%d_Tn%d_Kh%d_net%d_Gh%d_CI%d_sigma%0.2f_Bv%s_cost%s_mu%s_mig%s",N,TYPES,Tpos,Tneg,spar->kh,NETWORK,Gh,CI,spar->sigma,nparam[0],nparam[1],nparam[2],nparam[3]);
		#if (EVO==0)
		sprintf(nevo,"_dthmax%de%d_v0",mf,ef);
		#else
		sprintf(nevo,"_tlp");
		#endif
		#if (Tneg>0)
		sprintf(ntneg,"_CRnnA%d_CRnnB%0.1f_CRnpA%d_CRnpB%0.1f",(int)CRnn0,CRnn1,(int)CRnp0,CRnp1);
		#else
		sprintf(ntneg,"");
		#endif
	sprintf(meas->ftname_pars,"%s%s%s",ngeral,nevo,ntneg);
	free(ngeral);
	free(nevo);
	free(ntneg);
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
/*****************************************************
*   Populates Host and Microbial layer               *
*   Bacterial layer: types of bacteria               *
*   start with fixed fractions equal to 1/TYPES      *
******************************************************/
void initialStateFixedFrac(void){
       int i,j,nh=0,ne=0;
        double p_oc;

	p_oc=(double)H0/N;

	for(i=0; i<N; ++i){
		listh->vec[i]=0;
		inverselisth[i]=0;
	}

        for(i=0; i<N; ++i){
		if(FRANDOM<p_oc){
                        host[i]=1;
			listh->vec[nh]=i;
			inverselisth[i]=nh;
			++nh;
			spar->micr[i]=Bac0;
                        for(j=0; j<TYPES; ++j){
                                bac[i][j]=Bac0/TYPES;
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
/*****************************************************
*   Populates Host and Microbial layer               *
*   Bacterial layer: types of bacteria               *
*   are uniformly distributed                        *
*****************************************************/
void initialStateUniD(void){
       int i,j,nh=0,ne=0;
        double norm,p_oc;

	for(i=0; i<N; ++i){
		listh->vec[i]=0;
		inverselisth[i]=0;
	}

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
	
	for(i=0; i<N; ++i){
		listh->vec[i]=0;
		inverselisth[i]=0;
	}
                       
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
* (paper version, for TYPES=100)                 *
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
* (same as the paper version when TYPES=2,       *
* with 1 neutral and 1 positive)                 *
**************************************************/
void setInvestments(void){
        int i,neutral,ty_maior;
	double dty;

#if((Tpos+Tneg>TYPES)||(Tpos+Tneg<TYPES-1))
	printf("Tpos+Tneg should be equal to TYPES-1 (if you want to include a neutral type of bacteria) or TYPES (for no neutral type of bacteria).\n");
	exit(1);
#endif
#if(Tpos>=Tneg)
	ty_maior=Tpos;
#else 
	ty_maior=Tneg;
#endif
	dty=1./ty_maior;
	neutral=TYPES-Tpos-Tneg;/*=0 if there is no neutral type and =1 otherwise*/
	for(i=0; i<Tneg; ++i){
		spar->inv[i]=(double)(i-Tneg)/ty_maior;
		spar->s[i]=-1.;
	}
	for(i=Tneg; i<=Tpos+Tneg; ++i){
		spar->inv[i]=(double)(i-Tneg)/ty_maior + dty*(1.-(double)neutral);/*if there is a neutral type, its index is i=Tneg, and inv[Tneg]=0.. otherwise this is the index 
										   *if there is no neutral type, index i=Tneg is the index of the first positive type, and its investiment is inv[Tneg]=dty*/
		spar->s[i]=1.;
	}

        return;
}
/****************************************************************************
 *                     Set Initial Conditions                               *
 ***************************************************************************/
void setCI(TimeMeasures *meas){

	meas->Tnow=0.;
	meas->NTnow=0;
	/*setting initial state (alive hosts and bacteria abundances)*/
	#if (CI==0)//uniform distribution
		initialStateUniD();
	#elif(CI==1)//frequencies come from normal distribution
		initialStateNormD();
	#elif(CI==2)//single host
		initialStateSingleH(meas);
	#else
		initialStateFixedFrac();
	#endif
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
	setCI(meas);

	return;
}
