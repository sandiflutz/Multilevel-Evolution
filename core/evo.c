/* evo.c */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"globals.h"
#include"randgen_ufrgs.h"
#include"tools.h"
#include"evo.h"
#include"bac.h"
#include"measures.h"

/****************************************************************
*     calculates the accumulated investment in host @index      *
*     (related to host birth rate)                              *
*****************************************************************/
double calcAcumInvest(int index){
        int j;
	double cinv,fneg,kbac=spar->kbac;

        cinv=0.;

       	for(j=0; j<TYPES; ++j){
		cinv+=bac[index][j]*spar->inv[j]/kbac;
        }

        return cinv;
}
/********************************************************
*       Set host birth-death rates                      *
*********************************************************/
void setIndividualHostRates(Event *event){
        int i,idh,nh;
	int gh=spar->gh;
	int kh=spar->kh;
	double beta=Beta;
	double sb=Sb;
	double sd=Sd;
	double *w=NULL;

	nh=listh->usize;//number of hosts

	w=(double *)calloc(nh,sizeof(double));//vector for the accumulated investments of each live host

        event->sizeE=(2*nh);//# of possible host events (birth or death for each host)

	event->ratesE=(double *)calloc(event->sizeE,sizeof(double));

        //birth events
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		w[i]=calcAcumInvest(idh);
		event->ratesE[i]=(double)beta*(1.+sb*w[i])/gh;
	}

        //death events
	for(i=nh; i<(2*nh); ++i){
		event->ratesE[i]=beta*(1.-sd*w[i-nh])*nh/((double)kh*gh);
	}


	free(w);
        return;
}
/****************************************************************
*   For non well-mixed cases: 					*
*   each host event rate, associated to host @i, is 		*
*   substituted by the average host rate in the group centered	*
*   on @i. Probabilities are also multiplied by:		*
*   	-Repr. of i: rho_e[i]					*
*   	-Death: (1-rho_e[i])					*
*   	(rho_e[i]=fraction of empty sites in @i's group)	*
*****************************************************************/
void setGrRates(Event *event,Event *mevent){
	int i,k,idh,idviz,idlistv,nh,na;
	double *effective_rate = NULL;

	nh=listh->usize;

	effective_rate=(double *)calloc(2*nh,sizeof(double));

	if(spar->mh>0.){
		mevent->sizeE=nh;
		mevent->ratesE=(double *)calloc(mevent->sizeE,sizeof(double));
	}

	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		effective_rate[i]=event->ratesE[i];
		na=1;
		for(k=0; k<VIZ; ++k){
			idviz=neighbor[idh][k];
			if(host[idviz]==1){
				idlistv=inverselisth[idviz];
				effective_rate[i]+=event->ratesE[idlistv];
				++na;
			}
		}
		effective_rate[i]/=(double)na;
	}

	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		//event->ratesE[i]=effective_rate[i]*rho_e[idh];
		event->ratesE[i]=effective_rate[i]*ceil(rho_e[idh]);//if there is no empty sites on idh's neighborhood, its birth rate is 0
		if(spar->mh>0.){
			mevent->ratesE[i]=effective_rate[i]*(1.-rho_e[idh]);
		}
	}
	/*
	for(i=nh; i<2*nh; ++i){
		idh=listh->vec[i-nh];
		event->ratesE[i]*=(1.-rho_e[idh]);
	}
	*/
		

	free(effective_rate);
	return;
}
/****************************************************************
*       set cumulative rates for host events                    *
*****************************************************************/
void setCumulativeRates(Event *event){
	
	event->cprobE=(double *)calloc(event->sizeE,sizeof(double));
	cumulProb(event->sizeE,event->ratesE,event->cprobE);

	return;
}
/****************************************************************
*     Set microbial frequencies for the offspring of host @idp 	*
*     Bacteria types and their frequencies are randomly        	*
*     selected using a normal distribution around the bacteria	*
*     frequencies on the parent.                               	*
*****************************************************************/
void setMicrKidsNorm(int idp, int idk){
        int i,id,ns,p;
        double nr,fk,fp,norm,*cprob=NULL;
	double std=spar->sigma;
	double *gaussample;

	gaussample = (double *)calloc(SAMPLE,sizeof(double));

	cprob=(double *)calloc(TYPES,sizeof(double));//cumulative probability for the bacteria types
	cprob[0]=bac[idp][0]/spar->micr[idp];
		
	for(i=1; i<TYPES; ++i){
		cprob[i]=cprob[i-1]+bac[idp][i]/spar->micr[idp];
	}
		
	norm=0.;
	ns=0;
	while(ns<BSAMPLES){
		nr=FRANDOM*cprob[TYPES-1];
		p=selectEventCP(nr,cprob,TYPES);
		fp=bac[idp][p]/spar->micr[idp];
		normalRandSample(fp,std,0.,1.,SAMPLE,gaussample);
		id=(int)(FRANDOM*SAMPLE);
		fk=gaussample[id];
		bac[idk][p]+=fk;
		norm+=fk;
		++ns;
	}
	
	free(cprob);
	for(i=0; i<TYPES; ++i){
		bac[idk][i]*=Bacv/norm;//normalizing bac so the sum of microbes in host idk is bv
	}

	spar->micr[idk]=Bacv;

	free(gaussample);
        return;
}
/****************************************************************
*     Set microbial frequencies for the offspring of host @idp 	*
*     Bacteria types and their frequencies are randomly        	*
*     selected using a normal distribution around the bacteria	*
*     frequencies on the parent. 				*
*     -In this version, the amount of 				*
*     bacteria passed to the children may disappear from	*
*     the parent. 						*
*     -The amount of bacteria passed is a fixed frequency Fp 	*
*     of the amount of bacteria in the parent.			*
*     -Types of bacteria, in the parent, that have a frequency	*
*     of fp<eps are not included the kids microbiome		*
*****************************************************************/
void setMicrKidsNormPass(int idp, int idk){
	int id,j;
	double fpj,fkj,nk,norm;
	double np=spar->micr[idp];
	double std=spar->sigma;
	double *gaussample;

	gaussample = (double *)calloc(SAMPLE,sizeof(double));
	
	#if (WLMicr==0)
	nk=Bacv;
	#else//in this case part of the parent's microbes are give away to their children
	double fv=Fvert;
	nk=np*fv;
	#endif
	norm=0.;
	for(j=0; j<TYPES; ++j){
		fpj=bac[idp][j]/np;
		if(fpj>Fmin){
			normalRandSample(fpj,std,0.,1.,SAMPLE,gaussample);
			id=(int)(FRANDOM*SAMPLE);
			fkj=gaussample[id];
			norm+=fkj;
			bac[idk][j]=fkj;
		}
	}
	

	#if (WLMicr==0)
	for(j=0; j<TYPES; ++j){
		bac[idk][j]*=nk/norm;
	}
	#else//part of the parent's microbes are give away to their children
	spar->micr[idp]=0.;
	for(j=0; j<TYPES; ++j){
		bac[idk][j]*=nk/norm;
		bac[idp][j]-=bac[idk][j];
		spar->micr[idp]+=bac[idp][j];
	}
	#endif

	spar->micr[idk]=nk;

	free(gaussample);
	return;
}
/***************************************************************
*     Set microbial frequencies for the offspring of hosr @idp *
*     Bacteria types and their frequencies are randomly        *
*     selected using a poisson distribution for the number of  *
*     bacteria types in the offspring                          *
****************************************************************/
void setMicrKidsPoiss(int idp, int idk){
        int i,ns,p,sample;
        double nr,mean,*cprob=NULL;

        cprob=(double *)calloc(TYPES,sizeof(double));//cumulative probability for the bacteria types
	cprob[0]=bac[idp][0]/spar->micr[idp];
        for(i=1; i<TYPES; ++i){
                cprob[i]=cprob[i-1]+bac[idp][i]/spar->micr[idp];
        }

	mean=1.;
	ns=poissonRandKnuth(mean);
	sample=(int)fmax(1.,(double)ns);
        ns=0;
        while(ns<sample){
		nr=FRANDOM*cprob[TYPES-1];
		p=selectEventCP(nr,cprob,TYPES);
                bac[idk][p]+=(double)Bacv/sample;
		++ns;
        }
	spar->micr[idk]=Bacv;

        free(cprob);
        return;
}
/******************************************************
 *      birth of a new host                           *
 *****************************************************/
void hostBirth(int idp, int idk){

        host[idk]=1;
	#if (TV==0)//vertical transmission using a normal dist.
	setMicrKidsNorm(idp,idk);//set new host microbiome
	#elif (TV==1)//vertical transmission using a normal dist. with parent host loosing a fraction of their microbiome
	setMicrKidsNormPass(idp,idk);
	#else
	setMicrKidsPoiss(idp,idk);
	#endif


	/****measuraments****/
	#ifdef DIFBACOMPxT
		#if (OFFCOMP==0)
		storeBacDiffComp(idp,idk);
		#else
		offcomp->vecf[offcomp->usizef]=calcAcumInvest(idk)/spar->micr[idk];
		++offcomp->usizef;
		#endif
	#endif

        return;
}
/******************************************************
*       death of a host                               *
*******************************************************/
void hostDeath(int idh){
        int i;

        //host[idh]=0;
        host[idh]=-1;
        for(i=0; i<TYPES; ++i){
                bac[idh][i]=0.;
        }
	spar->micr[idh]=0.;

        return;
}
/****************************************
*	moviment of hosts:		*
*	2 neighbors exchange places	*
*****************************************/
void hostMoviment(int id1,int id2){
	int i,htmp;
	double bac_tmp;

	//microbiome
	for(i=0; i<TYPES; ++i){
		bac_tmp=bac[id1][i];
		bac[id1][i]=bac[id2][i];
		bac[id2][i]=bac_tmp;
	}
	//host positions
	htmp=host[id1];
	host[id1]=host[id2];
	host[id2]=htmp;

        return;
}
/************************************************
*         returns gillespie's time increment    *
************************************************/
double gillespieTime(double sumprob){
        double nr,passive_time;

	nr=FRANDOM;

        passive_time=-log(nr)/sumprob;

        return passive_time;
}
/**************************************************
*      Adjust host time step: the probability of  *
*      2 host events in a host timestep is <0.01  *
***************************************************/
double adjustTimeStep(double maxprob){
	int i,ok;
	double dt,*p0,*p1,*p2,max_p2=0.01;
	
	p0=(double *)calloc(DTVSIZE,sizeof(double));
	p1=(double *)calloc(DTVSIZE,sizeof(double));
	p2=(double *)calloc(DTVSIZE,sizeof(double));
	
	for(i=0; i<DTVSIZE; ++i){
		/*dtVec[] is a vector with a range of @DTVSIZE dt possible values in a logarithmic scale*/
		p0[i]=exp(-maxprob*dtVec[i]);
		p1[i]=maxprob*dtVec[i]*p0[i];
		p2[i]=1.-p0[i]-p1[i];
	}
	
	i=DTVSIZE-1;
	ok=0;
	dt=dtVec[0];
	do{
		if(p2[i]<max_p2){
			dt=dtVec[i];
			ok=1;		
		}

		--i;
	}while(i>0&&ok==0);

	free(p0);
	free(p1);
	free(p2);
	return dt;
}
/****************************************************************
 *	Returns the number of host time steps per bacteria 	*
 *  	time step						*
 ****************************************************************/ 
int hostNTSPerBacNTS(Event *event){
	int dnumsteps;
	double dt;

	dt=adjustTimeStep(event->cprobE[event->sizeE-1]);
	dnumsteps=ceil(Dt_ref/dt);
	dt=Dt_ref/(double)dnumsteps;
	stime->dth=dt;

	return dnumsteps;
}
/****************************************************************
*	Update group vacancy frequency centered on site @idh	*
*	due to a change in its state:				*
*	if site @idh is occupied, means there was a decrease	*
*	in the number of empty sites in the group, otherwise,	*
*	there was an increase.					*
*****************************************************************/
void updateEmptySpaceGrFreq(int idh){
	int k,idviz,sign;
	

	if(host[idh]==1){//if site is occupied, fraction of empty sites decreases 
		sign=-1;
	}else{//otherwise, fraction of empty sites in the group increases
		sign=1;
	}

	for(k=0; k<VIZ; ++k){
		idviz=neighbor[idh][k];
		rho_e[idviz]+=(double)sign*1./VIZ;
	}

	return;
}
/********************************************************
*	Choose site for migration:			*
*	states of sites @idm and one chosen are 	*
*	exchanged. Choice randomly chooses a 		*
*	neighboring site, according to probabilities 	*
*	that depend on vacancy (VIZ=#of neighbors):	*
*		-prob[i]=1/VIZ, if neighbor is empty,	*
*		-prob[i]=rho_e[i]/VIZ, otherwise	* 
*********************************************************/
int chooseMigSite(int idm){
	int i,k,id,idv;
	int *which_host = NULL;
	double *prob = NULL;    

	prob=(double *)calloc(VIZ,sizeof(double));
	which_host=(int *)calloc(VIZ,sizeof(int));

	for(i=0; i<VIZ; ++i){
		k=neighbor[idm][i];
		which_host[i]=k;
		if(host[k]==0){
			prob[i]=1./(double)VIZ;
		}else{
			prob[i]=rho_e[k]/(double)VIZ;
		}
	}
	id=selectEvent(FRANDOM,prob,VIZ);
	idv=which_host[id];


	free(prob);
	free(which_host);

	return idv;
}
/****************************************************************
*	Kill hosts that have a microbiome extremely low		*
*****************************************************************/
void killHostWithoutMicr(void){
	int i,idh,nh;

	nh=listh->usize;
	for(i=0; i<nh; ++i){
		nh=listh->usize;
		idh=listh->vec[i];
		if((spar->micr[idh]<EPS)&&(host[idh]==1)){
			host[idh]=0;
			spar->micr[idh]=0.;
			exchange(inverselisth,idh,listh->vec[nh-1]);
			exchange(listh->vec,i,nh-1);
			--listh->usize;

		}
	}

	return;
}
/************************************************
*	Dynamics for host migration events 	*
*	that happen in a Dt_ref (=microbial	*
*	time step)				*
*************************************************/
void hostMigrationDynamics(int dnumsteps,Event *mevent){
	int i,idh,idlist,idk,idlistk,nh,nm;
	double nr,pm;
	double mh=spar->mh;
	int *listm = NULL;
	
	nh=listh->usize;
	listm=(int *)calloc(dnumsteps,sizeof(int));
	nm=0;

	for(i=0; i<dnumsteps; ++i){
		nr=FRANDOM;
		pm=mevent->cprobE[nh-1]*mh*stime->dth;
		if(nr<pm){
			idlist=selectEventCP(nr,mevent->cprobE,nh);//selects a host for future migration
			listm[nm]=idlist;//storing its listh's id to keep track of its state and not the physical position
			++nm;
		}
	}
	                
	for(i=0; i<nm; ++i){//host migration list
		idlist=listm[i];
		idh=listh->vec[idlist];
		idk=chooseMigSite(idh);
		idlistk=inverselisth[idk];

		hostMoviment(idh,idk);//idh and idviz change places
		
		exchange(inverselisth,idh,idk);
		exchange(listh->vec,idlist,idlistk);
		if(host[idh]==0){//if idk was empty before the change, idh is now empty both of their groups suffer a change in the fraction of empty sites
			updateEmptySpaceGrFreq(idh);
			updateEmptySpaceGrFreq(idk);
		}

	}

	free(listm);
	return;
}
/****************************************************************
*	Host Layer Evolution (for a time interval=Dt_ref) for	*
*	the complete graph version				*
*****************************************************************/
void evolveHostCG(int dnumsteps,Event *event){
	int i,idh,idlist,ide,idk,nh,nb,nd,ne,whichE;
	double nr;
	int *listb = NULL;
	int *listd = NULL;
	DynList empty_viz;
	
	nh=listh->usize;
	
	empty_viz.size=N-nh;
	empty_viz.usize=0;
	empty_viz.vec = (int *)calloc(empty_viz.size,sizeof(int));
	for(i=0; i<empty_viz.size; ++i){
		ide=listh->vec[i+nh];
		empty_viz.vec[i]=listh->vec[i+nh];
		++empty_viz.usize;
	}
                
	listb=(int *)calloc(dnumsteps,sizeof(int));
	listd=(int *)calloc(dnumsteps,sizeof(int));
	nb=0;
	nd=0;

	for(i=0; i<dnumsteps; ++i){
		nr=FRANDOM;
		if(nr<event->cprobE[event->sizeE-1]){
			whichE=selectEventCP(nr,event->cprobE,event->sizeE);
			idh=listh->vec[whichE%nh];
			if(host[idh]==1){//if chosen host is alive and the system has more than 1 host
				switch(whichE/nh){
					case 0: 
						ne=empty_viz.usize;
						if(ne>0){
							ide=randNeighborID(idh,empty_viz.vec,0,ne,empty_viz.size);//randomly chooses an index of an empty site, stored on empty_viz.vec
							idk=empty_viz.vec[ide];//host id of the empty site is ide'th element of empty_viz.vec
							hostBirth(idh,idk);//site @idk receives offspring of @idh
							listb[nb]=idk;//@idk is stored on the list of newhosts (during this time interval of Dt_ref)
							++nb;//number of new hosts
							exchange(empty_viz.vec,ide,ne-1);
							--empty_viz.usize;//decreasing the number of available empty sites
							//measures	
							#ifdef GENTIME
							timeb[idk]=stime->Tnow+Dt_ref;//time of birth of host @idk
							if(timeb[idh]>=0.){
								meas->timegh+=(stime->Tnow+Dt_ref-timeb[idh]);//host @idh first reproduction time is the current time minos the time of its birth
								timeb[idh]=-1.;//just the first reproduction of an individual counts
								++meas->ngh;
							}
							#endif
						}
						break;
					case 1: 
						hostDeath(idh);
						listd[nd]=idh;//list of sites that received offspring during the Dt_ref time interval
						++nd;
						#ifdef GENTIME
						timeb[idh]=-1.;
						#endif
						break;
				}
			}
		}
	}
	                
	/*updating hosts dynamic list*/
	for(i=0; i<nb; ++i){
		nh=listh->usize;
		idh=listb[i];
		idlist=inverselisth[idh];
		exchange(inverselisth,listh->vec[nh],idh);
		exchange(listh->vec,nh,idlist);
		++listh->usize;
	}
	for(i=0; i<nd; ++i){
		nh=listh->usize;
		idh=listd[i];
		idlist=inverselisth[idh];
		exchange(inverselisth,idh,listh->vec[nh-1]);
		exchange(listh->vec,idlist,nh-1);
		--listh->usize;
		host[idh]=0;
	}
	
	/**************/
	#ifdef NUMHEVENTSxT
	meas->numb=nb;
	meas->numd=nd;
	#endif

	free(empty_viz.vec);
	free(listb);
	free(listd);
	return;
}
/****************************************************************
*	Host Layer Evolution (for a time interval=Dt_ref) for	*
*	the lattice version					*
*****************************************************************/
void evolveHostSL(int dnumsteps,Event *event){
	int i,idh,idlist,idk,ne,ide,nh,nb,nd,whichE;
	double nr;
	DynList empty_viz;
	double currentime=stime->Tnow+Dt_ref;
	
	nh=listh->usize;
	empty_viz.size=VIZ;
	empty_viz.usize=0;
	empty_viz.vec = (int *)calloc(empty_viz.size,sizeof(int));
                
	int *listb=(int *)calloc(dnumsteps,sizeof(int));
	int *listd=(int *)calloc(dnumsteps,sizeof(int));
	nb=0;
	nd=0;

	for(i=0; i<dnumsteps; ++i){
		nr=FRANDOM;
		if(nr<event->cprobE[event->sizeE-1]){
			whichE=selectEventCP(nr,event->cprobE,event->sizeE);
			idh=listh->vec[whichE%nh];
			if(host[idh]==1){//if chosen host is alive and the system has more than 1 host
				switch((int)(whichE/nh)){
					case 0: //birth 
						searchEmptyNeighbors(0,idh,host,neighbor,&empty_viz);
						ne=empty_viz.usize;
						if(ne>0){
							ide=randNeighborID(idh,empty_viz.vec,0,ne,empty_viz.size);//randomly chooses the index of an empty site, stored on empty_viz.vec
							idk=empty_viz.vec[ide];//host id of the empty site, idk, is the ide-th element of empty_viz.vec
							listb[nb]=idk;//list of sites that receive offspring during the Dt_ref time interval
							++nb;
							hostBirth(idh,idk);
							//measures	
							#ifdef GENTIME
							meas->timeR[idk][nR[idk]] = currentime;//time of birth of host @idk (this event is happening in Tnow<time<Tnow+Dt_ref)
							++nR[idk];
							
							meas->timeR[idh][nR[idh]] = currentime - meas->timeR[idh][nR[idh]-1];
							if(nR[idh]==1){//first reproduction
								timegh+=meas->timeR[idh][nR[idh]];
								++ngh;
							++nR[idh];
							}
							#endif
						}
						break;
					case 1: 
						listd[nd]=idh;
						++nd;
						hostDeath(idh);
						#ifdef GENTIME
						for(j=0; j<nR[idh]; ++j){
							meas->timeR[idh][j]=0.;
						}
						nR[idh]=0;
						#endif
						break;
				}
			}
		}
	}
	                
	/*updating hosts dynamic list*/
	for(i=0; i<nb; ++i){//reproduction list
		nh=listh->usize;
		idh=listb[i];
		idlist=inverselisth[idh];
		//updating group vacancy fraction (for groups around the newborn)
		updateEmptySpaceGrFreq(idh);
		//update lists related to positions of alive hosts
		exchange(inverselisth,listh->vec[nh],idh);
		exchange(listh->vec,nh,idlist);
		++listh->usize;
	}
	for(i=0; i<nd; ++i){//death list
		nh=listh->usize;
		idh=listd[i];
		idlist=inverselisth[idh];
		//updating group vacancy fraction (for groups around the dead host)
		updateEmptySpaceGrFreq(idh);
		//update lists related to positions of alive hosts
		exchange(inverselisth,idh,listh->vec[nh-1]);
		exchange(listh->vec,idlist,nh-1);
		--listh->usize;
		//official death is declared
		host[idh]=0;
	}
	
	/**************/
	#ifdef NUMHEVENTSxT
	meas->numb=nb;
	meas->numd=nd;
	#endif

	free(empty_viz.vec);
	free(listb);
	free(listd);
	return;
}
/****************************************
*          general time loop            *
*****************************************/
void callSysDynamics(Event *event, Event *mevent){
        int i,numsteps,nh,dnumsteps;

        numsteps=0;
        nh=listh->usize;
        dnumsteps=1;
	
        while((stime->Tnow<=stime->Tf)&&(nh>0)){
                #ifdef TMEAS
                timeMeasures();
                        #ifdef NUMHEVENTSxT
                        meas->numb=0;
                        meas->numd=0;
                        #endif
                #endif
                #ifdef STEADY_STATE_MEAS
                        if((stime->Tnow<stime->saveT+stime->timewindow)&&(stime->Tnow>=stime->saveT)){
                                avinv->vecf[avinv->usizef]=calcAverInv();
                                ++avinv->usizef;
				printf("rh=%f time=%f avinv[%d]=%f\n",(double)spar->kh/N,stime->Tnow,avinv->usizef-1,avinv->vecf[avinv->usizef-1]);
				stime->saveT+=10.;
                        }
                #endif
		
		setIndividualHostRates(event);
                #if (EVO==0)//complete graph with adjustable host dt
		setCumulativeRates(event);
		dnumsteps=hostNTSPerBacNTS(event);
		for(i=0; i<event->sizeE; ++i){
			event->cprobE[i]*=stime->dth;
		}
                evolveHostCG(dnumsteps,event);
		freeVecsEvent(event);

                #elif (EVO==1)//square lattice with adjustable host dt
		setGrRates(event,mevent);//group rates for birth, death and migration events
		setCumulativeRates(event);
		dnumsteps=hostNTSPerBacNTS(event);
		if(spar->mh>0.){
			freeVecsEvent(event);
			setCumulativeRates(mevent);//cumulative rates for migration events
			hostMigrationDynamics(dnumsteps,mevent);
			freeVecsEvent(mevent);
			//update host rates and probs (since maxprox prob. doesn't change, time substep is the same)
			setIndividualHostRates(event);
			setGrRates(event,mevent);
			setCumulativeRates(event);
		}
                for(i=0; i<event->sizeE; ++i){
			event->cprobE[i]*=stime->dth;
                	
		}
                evolveHostSL(dnumsteps,event);
			
		freeVecsEvent(event);
                #endif

                nh=listh->usize;

                evoBac(Dt_ref,stime->Tnow);
	//	killHostWithoutMicr();

                stime->Tnow+=Dt_ref;
                numsteps+=dnumsteps;
        }

        return;
}
/****************************************
*          1 host time loop             *
*****************************************/
void callSysDynamics1H(Event *event){
        int numsteps;

        numsteps=0;
        while(stime->Tnow<=stime->Tf){
                #ifdef TMEAS
                timeMeasures();
                #endif
                evoBac(Dt_ref,stime->Tnow);
                stime->Tnow+=Dt_ref;
                ++numsteps;
        }

        return;
}
/****************************************
*  Free allocated memory for vectors	*
*  that are part of an Event struct	*
*****************************************/ 
void freeVecsEvent(Event *event){

	free(event->ratesE);
	free(event->cprobE);
	event->ratesE=NULL;
	event->cprobE=NULL;

	return;
}	
