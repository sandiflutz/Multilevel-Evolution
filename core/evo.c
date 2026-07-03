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
*       Set host birth, death and migration rates       *
*********************************************************/
void setIndividualHostRates(Event *event,Event *mevent){
        int i,idh,nh;
	int gh=spar->gh;
	int kh=spar->kh;
	double *w,averagew,ntot,aver_rhoe;
	double beta=Beta;
	double sb=Sb;
	double sd=Sd;

	nh=listh->usize;//number of hosts

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

	w=(double *)calloc(nh,sizeof(double));
        //birth and death events
	averagew=0.;
	ntot=0.;
	aver_rhoe=0.;
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		w[i]=calcAcumInvest(idh);
		averagew+=w[i];
		ntot+=spar->micr[idh];
		aver_rhoe+=rho_e[idh];
		event->ratesE[i]=(double)beta*(1.+sb*w[i])/gh;
		event->ratesE[i+nh]=beta*(1.-sd*w[i])*nh/((double)kh*gh);//death events
	}
	aver_rhoe/=(double)nh;
	averagew/=ntot;
#if (NETWORK!=0)
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
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
	#if (MIGRATION_TYPE==0)//homogeneous host migration rate (no local component)
		mevent->ratesE[i]=spar->mh*beta*(1.+sb*averagew)*(1.-aver_rhoe)/gh;
	#elif (MIGRATION_TYPE==1)//local component comes from the host's neighborhood dilution
		mevent->ratesE[i]=spar->mh*(1.+sb*averagew)*(1.-rho_e[idh])/gh;
	#elif (MIGRATION_TYPE==2)//local component comes from the host's investment
		mevent->ratesE[i]=spar->mh*(1.+sb*w[i])*(1.-aver_rhoe)/gh;
	#else//local component comes from the host's neighborhood dilution its investment
		mevent->ratesE[i]=spar->mh*(1.+sb*w[i])*(1.-rho_e[idh])/gh;
	#endif
		if(rho_e[idh]==0.){
			event->ratesE[i]=0.;
		}
	}
#endif

	free(w);
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
*	2 states in the network 	*
*	exchange places			*
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
/****************************************************************************************
*	Choose a site to migrate within a specific connection distance			*
*****************************************************************************************/ 
int chooseMigSiteConDis(int id, int rf){
	int i,i0,k,idv,idlist,r,nlviz,viz;
	double sum_prob;
	DynList listviz;
	int *ilistviz;
	int *which_host = NULL;
	double *prob = NULL;    

	listviz.size=N;
	listviz.usize=0;
	listviz.vec=(int *)malloc(sizeof(int)*N);
	ilistviz=(int *)malloc(sizeof(int)*N);
	for(i=0; i<N; ++i){
		ilistviz[i]=-1;
	}

	r=0;
	nlviz=0;
	do{
		if(nlviz==0){
			for(i=0; i<con[id]; ++i){
				viz=neighbor[id][i];
				if(ilistviz[viz]==-1){
					listviz.vec[nlviz]=viz;
					ilistviz[viz]=nlviz;
					++nlviz;
				}
			}
			i0=0;
			++r;
			listviz.usize=nlviz;
		}else if(r<rf){
			for(i=i0; i<listviz.usize; ++i){
				for(k=0; k<con[i]; ++k){
					viz=neighbor[i][k];
					 if(ilistviz[viz]==-1){
						 listviz.vec[nlviz]=viz;
						 ilistviz[viz]=nlviz;
						 ++nlviz;
					 }
				}
			}
			++r;
			i0=listviz.usize;
			listviz.usize=nlviz;
		}
	}while(r<rf);


	which_host=(int *)malloc(sizeof(int)*nlviz);
	prob=(double *)malloc(sizeof(double)*nlviz);

	sum_prob=0.;
	for(i=0; i<nlviz; ++i){
		idv=listviz.vec[i];
		which_host[i]=idv;
		prob[i]=0.;
		if(host[idv]==0){
			prob[i]=1.;
		}else if(rho_e[idv]>0.){
			prob[i]=rho_e[idv];
		}
		sum_prob+=prob[i];
	}
	/*selecting a site based on prob[]*/
	if(sum_prob>0.){
		for(i=0; i<nlviz; ++i){
			prob[i]/=sum_prob;
		}
		idlist=selectEvent(FRANDOM,prob,nlviz);
	}else{//if the region is fully crowded, choose randomly with uniform dist.
		idlist=(int)(FRANDOM*nlviz);
	}

	
	free(ilistviz);
	free(listviz.vec);
	free(prob);
	free(which_host);
	return idlist;
}
/****************************************************************************************
* 		store site network id's for host migration and				*
* 		the related probability of being chosen 				*
* 		depending on dilution.							*
* 		Sites included are only the one at a given				*
* 		distance								*	
*****************************************************************************************/
double findMigSiteswithEuclRmig(int id0,int r,int i0,int imax,double da,int *which_host,double *prob){
	int i,idf;
	int xr,yr,x0,y0,xf,yf;
	double angle,anglerad,sum_prob;

	y0=(int)(id0/L);
	x0=id0-y0*L;

	sum_prob=0.;
	angle=0.;
	for(i=i0; i<imax; ++i){
		anglerad=((double)Pi)*(angle/180.);

		xr=(int)((double)r*cos(anglerad));
		yr=(int)((double)r*sin(anglerad));
		
		/*boundary conditions*/
		if(xr<0){//left
			xf=(x0+xr+L)%L;
		}else{//right
			xf=(x0+xr)%L;
		}
		
		if(yr<0){//up
			yf=(y0+yr+L)%L;
		}else{//down
			yf=(y0+yr)%L;
		}

		/*recovering network id and calculation its prob. of being chosen based on dilution*/
		idf=yf*L+xf;
		
		which_host[i]=idf;
		prob[i]=0.;
		if(host[idf]==0){
			prob[i]=1.;
		}else if(rho_e[idf]>0.){
			prob[i]=rho_e[idf];
		}
		sum_prob+=prob[i];

		angle+=da;
	}

	return sum_prob;
}
/****************************************************************
*	Choose site for migration within a distance r:		*
*	states of sites @idm and one chosen are exchanged. 	*
*	The probability of being chosen given that a site 	*
*	is within a distance r depends of on vacancy:		*
*		-prob[i]=1/sumprob, if neighbor is empty,	*
*		-prob[i]=rho_e[i]/sumprob, if site i not 	*
*			empty and rho_e[i]>0 for at least	* 
*			one of the sites being consider		*
*		-prob[i]=1/num_sites, for all sites being 	*
*		consider otherwise				*
*****************************************************************/
int chooseMigSite(int idm,int rmig){
	int i,i0,imax,idlist,id_site,nmig;
	int r,r0;
	double sum_prob,da;
	int *which_host = NULL;
	double *prob = NULL;    

	r0=1;
	nmig=0;
	for(r=1; r<=rmig; ++r){
		nmig+=VIZ*r;//number of angles to calculate the position of the sites that can be choose fo the jump
	}
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
	
	sum_prob=0.;
	i0=0;
	for(r=r0; r<=rmig; ++r){
		imax=r*VIZ+i0;
		da=360./((double)r*VIZ);
		sum_prob+=findMigSiteswithEuclRmig(idm,r,i0,imax,da,which_host,prob);
		i0=imax;
	}

	/*selecting a site based on prob[]*/
	if(sum_prob>0.){
		for(i=0; i<nmig; ++i){
			prob[i]/=sum_prob;
		}
		idlist=selectEvent(FRANDOM,prob,nmig);
	}else{//if the region is fully crowded, choose randomly with uniform dist.
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
*	Long-range random host migration:			*
*	randomly choose a site for migration under the  	*
*	the following constraints:				*
*		-its distance from the focus site (network 	*
*		id @idm) has to be equal or larger than a 	*
*		minimum rmin					*
*		-the fraction of empty sites in it 		*
*		neighborhood has to be larger than 0.		*
*****************************************************************/
int longRangeMigSite(int idm, int rmig,int rmin){
	int idf, ok,trials,maxtrials;
	int rx,ry,rf;
	int x0,y0,xf,yf,xmax,xmin,ymax,ymin;
	
	
	y0=(int)(idm/L);
	x0=idm-y0*L;
	trials=0;
	maxtrials=L;
	ok=0;
	do{
		idf=(int)(FRANDOM*N);
		yf=(int)(idf/L);
		xf=idf-yf*L;
		xmax=max(x0,xf);
		xmin=x0+xf-xmax;
		if((L-xmax+xmin)<(xmax-xmin)){//because of boundary conditions
			rx=L-xmax+xmin;
		}else{
			rx=xmax-xmin;
		}
		ymax=max(y0,yf);
		ymin=y0+yf-ymax;
		if((L-ymax+ymin)<(ymax-ymin)){//because of boundary conditions
			ry=L-ymax+ymin;
		}else{
			ry=ymax-ymin;
		}
		rf=(int)round(sqrt((rx*rx+ry*ry)));
		++trials;
		if((rf>=rmin)&&(idf!=idm)&&(rho_e[idf]>0.))ok=1;
		if(trials>maxtrials){//avoinding unending loop
			if(idf==idm){
				++maxtrials;
			}else{
				ok=1;
			}
		}
	}while(ok==0);

	return idf;
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
	int i,idh,idlist,idk,idlistk,nh,nm,ok;
	double nr,nrscaled,pm;
	int *listm = NULL;
	
	nh=listh->usize;
	listm=(int *)calloc(dnumsteps,sizeof(int));
	if(!listm){
		printf("It was not possible to allocate memory for vector listm.\n");
		exit(1);
	}
	nm=0;

	pm=mevent->cprobE[nh-1]*stime->dth;
	for(i=0; i<dnumsteps; ++i){
		nr=FRANDOM;
		if(nr<pm){
			#if (MIGRATION_TYPE==0)
			ok=0;
			do{
				idlist=(int)(FRANDOM*nh);
				if(host[listh->vec[idlist]]==1)ok=1;
			}while(ok==0);
			#else
			nrscaled=FRANDOM*mevent->cprobE[mevent->sizeE-1];
			idlist=selectEventCP(nrscaled,mevent->cprobE,mevent->sizeE);//selects a host for future migration
			#endif
			listm[nm]=idlist;//storing its listh's id to keep track of its state and not the physical position
			++nm;
		}
	}
	                
	for(i=0; i<nm; ++i){//host migration list
		idlist=listm[i];
		idh=listh->vec[idlist];//current position of host tagged as @idlist
	#if (NETWORK==1)
		#if (LONG_RANGE_MIG!=0)
		if(FRANDOM<spar->plr){
			idk=longRangeMigSite(idh,spar->rmigh,Rmin);
		}else{
			idk=chooseMigSite(idh,spar->rmigh);
		}
		#else
		idk=chooseMigSite(idh,spar->rmigh);
		#endif
	#elif (NETWORK==2)
		idk=chooseMigSiteConDis(idh,spar->rmigh);
	#endif
		idlistk=inverselisth[idk];

		hostMoviment(idh,idk);//idh and idk change places
		
		exchange(inverselisth,idh,idk);
		exchange(listh->vec,idlist,idlistk);
		if((host[idh]==0)&&(host[idk]!=0)){//if idk was empty before the change, idh is now empty. Both of their groups suffer a change in the fraction of empty sites (otherwise, nothing changes, 
						   //because idh has to be occupied before migration)
			updateEmptySpaceGrFreq(idh);
			updateEmptySpaceGrFreq(idk);
		}
	}

	#ifdef NUMHEVENTSxT
	meas->nummh+=nm;
	#endif
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
	double nr,nrscaled,pbd;
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
		pbd=event->cprobE[event->sizeE-1];
		if(nr<pbd){
			nrscaled=FRANDOM*event->cprobE[event->sizeE-1];
			whichE=selectEventCP(nrscaled,event->cprobE,event->sizeE);
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
	meas->numb+=nb;
	meas->numd+=nd;
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
	double nr,nrscaled,pbd;
	DynList empty_viz;
	
	nh=listh->usize;
	#if (NETWORK==1)
	empty_viz.size=VIZ;
	#elif (NETWORK==2)
	empty_viz.size=maxcon;
	#endif
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
		pbd=event->cprobE[event->sizeE-1];
		if(nr<pbd){
			nrscaled=FRANDOM*event->cprobE[event->sizeE-1];
			whichE=selectEventCP(nrscaled,event->cprobE,event->sizeE);
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
	meas->numb+=nb;
	meas->numd+=nd;
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
                       
	#ifdef NUMHEVENTSxT
	meas->numb=0;
	meas->numd=0;
	meas->nummh=0;
	#endif

        while((stime->Tnow<=stime->Tf)&&(finish==0)){
                #ifdef TMEAS
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
                        if((stime->Tnow>=stime->saveT-0.0001)&&(stime->Tnow<=stime->saveT+0.0001)){
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
				#ifdef COSTxPLRxAVINV
				printf("cost=%f plr=%f time=%f avinv[%d]=%f\n",spar->cost,spar->plr,stime->Tnow,avinv->usizef-1,avinv->vecf[avinv->usizef-1]);
				#endif
				#ifdef MHxPLRxAVINV
				printf("(Storing time) mh=%f plr=%f avinv[%d]=%f time=%f\n",spar->mh,spar->plr,avinv->usizef-1,avinv->vecf[avinv->usizef-1],stime->Tnow);
				#endif
				stime->saveT=stime->Tnow+stime->tinterval;
                        }else{
				#ifdef MHxPLRxAVINV
				printf("(Transient time) time=%f\n",stime->Tnow);
				#endif
			}
                #endif
		#ifdef MULTIPLE_COSTS_WxT
		timeMeasures();
		#endif
		
		setIndividualHostRates(event,mevent);
                #if (NETWORK==0)//well-mixed
		setCumulativeRates(event);
		dnumsteps=hostNTSPerBacNTS(event);
		for(i=0; i<event->sizeE; ++i){
			event->cprobE[i]*=stime->dth;
		}
                evolveHostCG(dnumsteps,event);
		freeVecsEvent(event);

                #else//square lattice or smallworld
		setCumulativeRates(event);
		dnumsteps=hostNTSPerBacNTS(event);
		if(spar->mh>0.){
			freeVecsEvent(event);
			setCumulativeRates(mevent);//cumulative rates for migration events
			hostMigrationDynamics(dnumsteps,mevent);
			freeVecsEvent(mevent);
			//update host rates and probs (since maxprox prob. doesn't change, time substep is the same)
			setIndividualHostRates(event,mevent);
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
