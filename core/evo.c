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

	if(host[index]==1){
		for(j=0; j<TYPES; ++j){
			cinv+=bac[index*TYPES+j]*spar->inv[j]/kbac;
		}
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
	if(!w){
		printf("It wasn't possible to allocate memory for vector w.\n");
		exit(1);
	}	

        event->sizeE=(2*nh);//# of possible host events (birth or death for each host)

	if(event->ratesE){
		free(event->ratesE);
		event->ratesE=NULL;
	}
	event->ratesE=(double *)calloc(event->sizeE,sizeof(double));
	if(!event->ratesE){
		printf("It wasn't possible to allocate memory for vector event->ratesE.\n");
		exit(1);
	}	

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


	if(w){
		free(w);
		w=NULL;
	}
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
	int i,idh,nh;

	nh=listh->usize;

	if(spar->mh>0.){
		mevent->sizeE=nh;
		if(mevent->ratesE){
			free(mevent->ratesE);
			mevent->ratesE=NULL;
		}
		mevent->ratesE=(double *)calloc(mevent->sizeE,sizeof(double));
		if(!mevent->ratesE){
			printf("It wasn't possible to allocate memory for mevent->ratesE.\n");
			exit(1);
		}
	}
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		if(spar->mh>0.){
			mevent->ratesE[i]=event->ratesE[i]*(1.-rho_e[idh]);
		}
		event->ratesE[i]*=ceil(rho_e[idh]);//if there is no empty sites on idh's neighborhood, its birth rate is 0
	}
	return;
}
/****************************************************************
*       set cumulative rates for host events                    *
*****************************************************************/
void setCumulativeRates(Event *event){
	
	if(event->cprobE){
		free(event->cprobE);
		event->cprobE=NULL;
	}
	event->cprobE=(double *)calloc(event->sizeE,sizeof(double));
	if(!event->cprobE){
		printf("It wasn't possible to allocate memory for mevent->cprobE.\n");
		exit(1);
	}

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
	double *gaussample=NULL;

	gaussample = (double *)calloc(SAMPLE,sizeof(double));
	if(!gaussample){
		printf("It wasn't possible to allocate memory for vector gaussample on setMicrKidsNorm().\n");
		exit(1);
	}

	cprob=(double *)calloc(TYPES,sizeof(double));//cumulative probability for the bacteria types
	if(!cprob){
		printf("It wasn't possible to allocate memory for vector cprob on setMicrKidsNorm().\n");
		exit(1);
	}
	cprob[0]=bac[idp*TYPES]/spar->micr[idp];
		
	for(i=1; i<TYPES; ++i){
		cprob[i]=cprob[i-1]+bac[idp*TYPES+i]/spar->micr[idp];
	}
		
	norm=0.;
	ns=0;
	while(ns<BSAMPLES){
		nr=FRANDOM*cprob[TYPES-1];
		p=selectEventCP(nr,cprob,TYPES);
		fp=bac[idp*TYPES+p]/spar->micr[idp];
		normalRandSample(fp,std,0.,1.,SAMPLE,gaussample);
		id=(int)(FRANDOM*SAMPLE);
		fk=gaussample[id];
		bac[idk*TYPES+p]+=fk;
		norm+=fk;
		++ns;
	}

	if(cprob){
		free(cprob);
		cprob=NULL;
	}

	for(i=0; i<TYPES; ++i){
		bac[idk*TYPES+i]*=Bacv/norm;//normalizing bac so the sum of microbes in host idk is bv
	}

	spar->micr[idk]=Bacv;

	if(gaussample){
		free(gaussample);
		gaussample=NULL;
	}
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
	if(!cprob){
		printf("It wasn't possible to allocate memory for vector cprob on setMicrKidsPois().\n");
		exit(1);
	}
	cprob[0]=bac[idp*TYPES]/spar->micr[idp];
        for(i=1; i<TYPES; ++i){
                cprob[i]=cprob[i-1]+bac[idp*TYPES+i]/spar->micr[idp];
        }

	mean=1.;
	ns=poissonRandKnuth(mean);
	sample=(int)fmax(1.,(double)ns);
        ns=0;
        while(ns<sample){
		nr=FRANDOM*cprob[TYPES-1];
		p=selectEventCP(nr,cprob,TYPES);
                bac[idk*TYPES+p]+=(double)Bacv/sample;
		++ns;
        }
	spar->micr[idk]=Bacv;

        if(cprob){
		free(cprob);
		cprob=NULL;
	}
        return;
}
/******************************************************
 *      birth of a new host                           *
 *****************************************************/
void hostBirth(int idp, int idk){

        host[idk]=1;
	#if (TV==0)//vertical transmission using a normal dist.
	setMicrKidsNorm(idp,idk);//set new host microbiome
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
                bac[idh*TYPES+i]=0.;
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
		bac_tmp=bac[id1*TYPES+i];
		bac[id1*TYPES+i]=bac[id2*TYPES+i];
		bac[id2*TYPES+i]=bac_tmp;
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
	if(!p0){
		printf("It wasn't possible to allocate memory for vector p0 on adjustTimeStep().\n");
		exit(1);
	}
	p1=(double *)calloc(DTVSIZE,sizeof(double));
	if(!p1){
		printf("It wasn't possible to allocate memory for vector p1 on adjustTimeStep().\n");
		exit(1);
	}
	p2=(double *)calloc(DTVSIZE,sizeof(double));
	if(!p2){
		printf("It wasn't possible to allocate memory for vector p2 on adjustTimeStep().\n");
		exit(1);
	}
	
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

        if(p0){
		free(p0);
		p0=NULL;
	}
        if(p1){
		free(p1);
		p1=NULL;
	}
        if(p2){
		free(p2);
		p2=NULL;
	}
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
int chooseMigSite(int idm,int rmig){
	int i,id,idlist,id_site;
	int nmig;
	int xr,yr,x0,y0,xf,yf;
	int *which_host = NULL;
	double sum_prob,angle,anglerad,da,cossine,sine;
	double *prob = NULL;    

	nmig=VIZ*rmig;

	prob=(double *)calloc(nmig,sizeof(double));
	if(!prob){
		printf("It wasn't possible to allocate memory for vector prob on chooseMigSite().\n");
		exit(1);
	}
	which_host=(int *)calloc(nmig,sizeof(int));
	if(!which_host){
		printf("It wasn't possible to allocate memory for vector which_host on chooseMigSite().\n");
		exit(1);
	}

	y0=(int)(idm/L);
	x0=idm-y0*L;

	da=360./((double)nmig);
	angle=0.;
	sum_prob=0.;
	for(i=0; i<nmig; ++i){
		anglerad=((double)Pi)*(angle/180.);
		sine=sin(anglerad);
		cossine=cos(anglerad);
		xf=x0;
		yf=y0;

		xr=(int)(round((double)rmig*cossine));
		if(xr<0){//left
			xf=(x0+xr+L)%L;
		}else if(xr>0){//right
			xf=(x0+xr)%L;
		}
		
		yr=(int)(round((double)rmig*sine));
		if(yr<0){//up
			yf=(y0+yr+L)%L;
		}else if(yr>0){//down
			yf=(y0+yr)%L;
		}

		id=yf*L+xf;
		
		which_host[i]=id;
		prob[i]=0.;
		if(host[id]==0){
			prob[i]=1.;
		}else if(rho_e[id]>0.){
			prob[i]=rho_e[id];
		}
		sum_prob+=prob[i];

		angle+=da;
	}
	
	if(sum_prob>0.){
		for(i=0; i<nmig; ++i){
			prob[i]/=sum_prob;
		}
		idlist=selectEvent(FRANDOM,prob,nmig);
	}else{
		idlist=(int)(FRANDOM*nmig);
	}
	id_site=which_host[idlist];

	if(prob){
		free(prob);
		prob=NULL;
	}
	if(which_host){
		free(which_host);
		which_host=NULL;
	}
	return id_site;
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
	int *listm = NULL;
	
	nh=listh->usize;
	listm=(int *)calloc(dnumsteps,sizeof(int));
	if(!listm){
		printf("It was not possible to allocate memory for vector listm.\n");
		exit(1);
	}
	nm=0;

	for(i=0; i<dnumsteps; ++i){
		nr=FRANDOM;
		pm=mevent->cprobE[nh-1]*spar->mh*stime->dth;
		if(nr<pm){
			idlist=selectEventCP(nr,mevent->cprobE,mevent->sizeE);//selects a host for future migration
			listm[nm]=idlist;//storing its listh's id to keep track of its state and not the physical position
			++nm;
		}
	}
	                
	for(i=0; i<nm; ++i){//host migration list
		idlist=listm[i];
		idh=listh->vec[idlist];
		idk=chooseMigSite(idh,spar->rmigh);
		idlistk=inverselisth[idk];

		hostMoviment(idh,idk);//idh and idviz change places
		
		exchange(inverselisth,idh,idk);
		exchange(listh->vec,idlist,idlistk);
		if(host[idh]==0){//if idk was empty before the change, idh is now empty both of their groups suffer a change in the fraction of empty sites
			updateEmptySpaceGrFreq(idh);
			updateEmptySpaceGrFreq(idk);
		}
	}

	if(listm){
		free(listm);
		listm=NULL;
	}
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
	
	nh=listh->usize;
	
	ne=N-nh;
	listb=(int *)calloc(dnumsteps,sizeof(int));
	if(!listb){
		printf("It was not possible to allocate memory for vector listb.\n");
		exit(1);
	}
	listd=(int *)calloc(dnumsteps,sizeof(int));
	if(!listd){
		printf("It was not possible to allocate memory for vector listd.\n");
		exit(1);
	}
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
						if(ne>0){
							do{
								ide=randNeighborID(idh,listh->vec,nh,N,N);//randomly chooses an index of an empty site, stored on empty_viz.vec
								idk=listh->vec[ide];//host id of the empty site is ide'th element of empty_viz.vec
							}while(host[idk]!=0);
							hostBirth(idh,idk);//site @idk receives offspring of @idh
							listb[nb]=idk;//@idk is stored on the list of newhosts (during this time interval of Dt_ref)
							++nb;//number of new hosts
							--ne;//decreasing the number of available empty sites
							//measures	
						}
						break;
					case 1: 
						hostDeath(idh);
						listd[nd]=idh;//list of sites that received offspring during the Dt_ref time interval
						++nd;
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

	if(listb){
		free(listb);
		listb=NULL;
	}
	if(listd){
		free(listd);
		listd=NULL;
	}
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
	
	nh=listh->usize;
	empty_viz.size=VIZ;
	empty_viz.usize=0;
	empty_viz.vec = (int *)calloc(empty_viz.size,sizeof(int));
	if(!empty_viz.vec){
		printf("It was not possible to allocate memory for vector empty_viz.vec.\n");
		exit(1);
	}
                
	int *listb=(int *)calloc(dnumsteps,sizeof(int));
	if(!listb){
		printf("It was not possible to allocate memory for vector listb.\n");
		exit(1);
	}
	int *listd=(int *)calloc(dnumsteps,sizeof(int));
	if(!listd){
		printf("It was not possible to allocate memory for vector listd.\n");
		exit(1);
	}
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
						}
						break;
					case 1: 
						listd[nd]=idh;
						++nd;
						hostDeath(idh);
						break;
				}
			}
		}
	}
	                
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

	if(empty_viz.vec){
		free(empty_viz.vec);
		empty_viz.vec=NULL;
	}
	if(listb){
		free(listb);
		listb=NULL;
	}
	if(listd){
		free(listd);
		listd=NULL;
	}
	return;
}
/****************************************
*          general time loop            *
*****************************************/
void callSysDynamics(Event *event, Event *mevent){
        int i,numsteps,nh,dnumsteps,finish=0;

        numsteps=0;
        nh=listh->usize;
        dnumsteps=1;

	sysmeas->averw=0.;
	sysmeas->averw2=0.;
	sysmeas->nw=0;
        while((stime->Tnow<=stime->Tf)&&(finish==0)){
                #ifdef TMEAS
                        #ifdef NUMHEVENTSxT
                        meas->numb=0;
                        meas->numd=0;
                        #endif
			#ifdef BESTCLUSTER_TIMES
			lb=(int *)calloc(listh->usize,sizeof(int));
			if(!lb){
				printf("It wasn't possible to allocate memory for vector lb on callSysDynamics()\n");
				exit(1);
			}
			memset(lb,-1,sizeof(int)*listh->usize);
			maxclw->whichLBF=0;
                        #endif
                	timeMeasures();
                #endif
                #ifdef STEADY_STATE_MEAS
                        if((stime->Tnow<stime->saveT+stime->timewindow)&&(stime->Tnow>=stime->saveT)){
                                avinv->vecf[avinv->usizef]=calcAverInv();
                                ++avinv->usizef;
				#ifdef AVINVxRH
				printf("rh=%f time=%f avinv[%d]=%f\n",(double)spar->kh/N,stime->Tnow,avinv->usizef-1,avinv->vecf[avinv->usizef-1]);
				#endif
				#ifdef AVINVxGH
				printf("gh=%d time=%f avinv[%d]=%f\n",spar->gh,stime->Tnow,avinv->usizef-1,avinv->vecf[avinv->usizef-1]);
				#endif
				#ifdef AVINVxCOST
				printf("cost=%f time=%f avinv[%d]=%f\n",spar->cost,stime->Tnow,avinv->usizef-1,avinv->vecf[avinv->usizef-1]);
				#endif
				stime->saveT+=10.;
                        }
                #endif
		
		setIndividualHostRates(event);
                #if (NETWORK==0)//complete graph with adjustable host dt
		setCumulativeRates(event);
		dnumsteps=hostNTSPerBacNTS(event);
		for(i=0; i<event->sizeE; ++i){
			event->cprobE[i]*=stime->dth;
		}
                evolveHostCG(dnumsteps,event);
		freeVecsEvent(event);

                #elif (NETWORK==1)//square lattice with adjustable host dt
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
		#ifdef BESTCLUSTER_TIMES
		if(lb){
			free(lb);
			lb=NULL;
		}
		#endif
                
		#endif

                nh=listh->usize;

                evoBac(Dt_ref,stime->Tnow);

                stime->Tnow+=Dt_ref;
                numsteps+=dnumsteps;
		if(nh==0)finish=1;
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

	if(event->ratesE){
		free(event->ratesE);
		event->ratesE=NULL;
	}
	if(event->cprobE){
		free(event->cprobE);
		event->cprobE=NULL;
	}

	return;
}	
