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
void allocateMemory(void){
        int i;
	int e0,ef,m0,mf;//variable related to the range of dtVec: [m0*10^e0;mf*10^ef]

	/*host network*/
        host=(int *)calloc(N,sizeof(int));
	memset(host,0,sizeof(int)*N);
	
	#if (NETWORK!=0)//not the well-mixed/complete graph case
	/*network array: e.g. neighbor[k][idh]=idh_viz (label of the k-th neighbor of host @idh is @idh_viz)*/
	neighbor=(int **)calloc(N,sizeof(int *));
	for(i=0; i<N; ++i){
                neighbor[i]=(int *)calloc(VIZ,sizeof(int));
        }
	/*vector that stores the frequency of empty sites in each group @i, centered on site @i */
	rho_e=(double *)calloc(N,sizeof(double));
        #endif
        
	/*microbial abundance matrix: bac[i:host index][j:type of bacteria index]*/
        bac=(double **)calloc(N,sizeof(double *));
	for(i=0; i<N; ++i){
        	bac[i]=(double *)calloc(TYPES,sizeof(double));
	}
	/*cost vector: proportional to the cost constant gamma*/
	#if (Tneg>0)
        costvec=(double **)calloc(N,sizeof(double *));
	for(i=0; i<N; ++i){
        	costvec[i]=(double *)calloc(TYPES,sizeof(double));
	}
	#endif
	
	/**/	
	inverselisth=(int *)calloc(N,sizeof(int));

	/*subtimestep vector, to adjust number of times steps for host events*/
	dtVec=(double *) calloc(DTVSIZE,sizeof(double));
	m0=1;
	e0=-7;
	mf=1;
	ef=-2;
	logSpacedVec(dtVec,e0,ef,m0,mf,DTVSIZE);/*filling dt vector @dtVec with values from m0*10^(e0) to mf*10^(ef) logarithmic spaced
					       *passing: 1)dt vector,2)e0,3)ef,4)m0,5)mf,6)size of the dt vector  
					       *range: dtVec[0]=1e-07 and dtVec[DTVSIZE-1]=1e-02 are values used in the paper
					       *another option: dtVec[DTVSIZE-1]=0.05 is the value of the microbial time step in the original paper
					       */
        
	/****structs****/
        
	//system parameters necessary for the equations of the microbial
        
	spar = malloc(sizeof(SysParams));
        if (!spar) { perror("malloc"); exit(1);}
	spar->kh=K_H;
	spar->gh=Gh;
	spar->mh=Mh;
	spar->kbac=K_bac;
	spar->mu=Mu;
	spar->cost=Gamma;
	spar->mig=Theta;
	spar->sigma=SIGMA;
	spar->micr=(double *)calloc(N,sizeof(double));
	spar->inv=(double *)calloc(TYPES,sizeof(double));

	//system times
	stime = malloc(sizeof(SysTimes));
        if (!stime) { perror("malloc"); exit(1);}
	stime->Tnow=0.;
	stime->Tf=TF;
	stime->tinterval=10.;
	stime->transtime=Ttrans;
	stime->timewindow=Twin;
	stime->dth=Dt_ref;

	//lists structs
        listh = malloc(sizeof(DynList));
        if (!listh) { perror("malloc"); exit(1);}
        listh->vec = (int *)calloc(N,sizeof(int));
        listh->size=N;
        listh->usize=0;

	//system measures
	sysmeas = malloc(sizeof(SysMeasures));
	sysmeas->averw=0.;
	sysmeas->averw2=0.;
	sysmeas->nw=0;

       	return;
}
/*****************************************************
*   Populates Host and Microbial layer               *
*   Bacterial layer: types of bacteria               *
*   are uniformly distributed                        *
*****************************************************/
void initialStateUniD(void){
       int i,j,k,idviz,nh=0,ne=0;
        double norm,p_oc;

	for(i=0; i<N; ++i){
		listh->vec[i]=0;
		inverselisth[i]=0;
	}

	p_oc=(double)H0/N;

	#if (NETWORK!=0)//not well-mixed
	for(i=0; i<N; ++i){
		rho_e[i]=0.;
	}
	#endif
        
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
			host[i]=0;
			listh->vec[N-1-ne]=i;//empty sites are stored at the end of the list of hosts;
			inverselisth[i]=N-1-ne;
			++ne;
			memset(bac[i],0.,sizeof(double)*TYPES);
			spar->micr[i]=0.;
			#if (NETWORK!=0)//not well-mixed
			for(k=0; k<VIZ; ++k){
				idviz=neighbor[i][k];
				rho_e[idviz]+=1./VIZ;
			}
			#endif
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
       int i,j,k,idviz,nh=0,ne=0;
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

	#if (NETWORK!=0)//not well-mixed
	for(i=0; i<N; ++i){
		rho_e[i]=0.;
	}
	#endif
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
                        host[i]=0;
			listh->vec[N-1-ne]=i;//empty sites are stored at the end of the list of hosts
			inverselisth[i]=N-1-ne;
			++ne;
			memset(bac[i],0.,sizeof(double)*TYPES);
			spar->micr[i]=0.;
			#if (NETWORK!=0)//not well-mixed
			for(k=0; k<VIZ; ++k){
				idviz=neighbor[i][k];
				rho_e[idviz]+=1./VIZ;
			}
			#endif
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
void initialStateSingleH(void){
       int i;
        double norm;

        
	for(i=0; i<N; ++i){
		memset(bac[i],0.,sizeof(double)*TYPES);
		spar->micr[i]=0.;
	}
	
	host[0]=1;
	spar->micr[0]=Bac0;
	norm=0.;
        for(i=0; i<TYPES; ++i){
                bac[0][i]=FRANDOM;
                norm+=bac[0][i];
        }
        for(i=0; i<TYPES; ++i){
                bac[0][i]*=Bac0/norm;
        }
        for(i=0; i<N; ++i){
                listh->vec[i]=i;
		inverselisth[i]=i;
        }

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
	}
	for(i=Tneg; i<=Tpos+Tneg; ++i){
		spar->inv[i]=(double)(i-Tneg)/ty_maior + dty*(1.-(double)neutral);/*if there is a neutral type, its index is i=Tneg, and inv[Tneg]=0.. otherwise this is the index 
										   *if there is no neutral type, index i=Tneg is the index of the first positive type, and its investiment is inv[Tneg]=dty*/
	}

        return;
}
/************************************************************************
* Set Cost Vector: 							*
* 	costvec[id_host][id_type]=cost*inv[id_type]*func    (Tneg>0)	*
*************************************************************************/
void setCostVec(void){
	int i,j;

	/*stting cost vector*/

	for(i=0; i<N; ++i){
		for(j=0; j<TYPES; ++j){
			costvec[i][j]=spar->cost*spar->inv[j];
		}
	}

	return;
}
/****************************************************************************
 *                     Set Initial Conditions                               *
 ***************************************************************************/
void setCI(void){

	stime->Tnow=0.;
	/*setting investment vector*/
	setInvestments();
	/*setting initial state (alive hosts and bacteria abundances)*/
	#if (CI==0)//uniform distribution
		initialStateUniD();
	#elif(CI==1)//frequencies come from normal distribution
		initialStateNormD();
	#else//single host
		initialStateSingleH();
	#endif
	#if (Tneg>0)
	setCostVec();
	#endif
	
	return;
}
/****************************************************************************
 *                     Build System                                         *
 ***************************************************************************/
void setSystem(void){
	
        (void)start_randomic(0);

	/*Allocating memory*/
	allocateMemory();

	/*setting hosts network*/	
	#if(NETWORK==1)
		squareLattice(neighbor,VIZ,L);//from tools
	#endif

	/*setting initial state (alive hosts and bacteria abundances)*/
	setCI();

	return;
}
