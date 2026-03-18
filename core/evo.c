/* evo.c */

#include<math.h>
#include"globals.h"
#include"randgen_ufrgs.h"
#include"tools.h"
#include"evo.h"

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
/****************************************************
* Calculates the host event rates: birth and death  *
*****************************************************/
void calcHostEvents(Event *event){
        int i,idh,nh;
	int gh=spar->gh;
	int kh=spar->kh;
	double beta=spar->beta;
	double sb=spar->sb;
	double sd=spar->sd;
	double *w=NULL;

	nh=listh->usize;//number of hosts

	w=(double *)calloc(nh,sizeof(double));//vector for the accumulated investments of each live host

        event->usizeE=(2*nh);//# of possible host events (birth or death for each host)
        memset(event->ratesE,0.,sizeof(double)*event->sizeE);

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

	//calculate cumulative rates (that are gonna be used in a bissection method to randomly select an event)
	cumulProb(event->usizeE,event->ratesE,event->cprobE);

	free(w);
        return;
}
/****************************************************************
*     Set microbial frequencies for the offspring of host @idp 	*
*     Bacteria types and their frequencies are randomly        	*
*     selected using a normal distribution around the bacteria	*
*     frequencies on the parent.                               	*
*****************************************************************/
void setMicrKidsNorm(int idp, int idk){
        int i,id,id1,id2,ns,p;
        double nr,fk,fp,norm,*cprob=NULL;
	double std=spar->sigma;
	double *gaussample;

	gaussample = (double *)calloc(GSAMPLE,sizeof(double));

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
		normalRandSample(fp,std,0.,1.,GSAMPLE,gaussample);
		id=(int)(FRANDOM*GSAMPLE);
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
	double fv=spar->fvert;
	double *gaussample;

	gaussample = (double *)calloc(GSAMPLE,sizeof(double));
	
	#if (WLMicr==0)
	nk=Bacv;
	#else//in this case part of the parent's microbes are give away to their children
	nk=np*fv;
	#endif
	norm=0.;
	for(j=0; j<TYPES; ++j){
		fpj=bac[idp][j]/np;
		if(fpj>fv){
			normalRandSample(fpj,std,0.,1.,GSAMPLE,gaussample);
			id=(int)(FRANDOM*GSAMPLE);
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

        ns=0;
	mean=10e03*Bacv;
	if(mean<=30.){
		sample=poissonRandKnuth(mean);
	}else{
		double ng=gaussRandNum(mean,mean);
		sample=round(ng);
		sample=(int)(fmax((double)sample,0.));
	}
        while(ns<sample){
		nr=FRANDOM*cprob[TYPES-1];
		p=selectEventCP(nr,cprob,TYPES);
                bac[idk][p]+=(double)Bacv/sample;
		++ns;
        }
	if(sample>0)spar->micr[idk]=Bacv;

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
/************************************************
*	After a host is chosen to 		*
*	reproduce, it has a probability		* 
*	to reproduce that depends on the	* 
*	number of empty sites available		*
*	to them. This routine returns 		*
*	this probability			*
*************************************************/
double birthFunc(int ne,int nemax,int whichfunc){
	double de,fb;
		
	if(nemax==0){
		printf("Maximum number of empty sites, passed to birthFunc(), cannot be 0.");
		exit(1);
	}else{
		de=(double)ne/nemax;
	}

	if(whichfunc==0){//no host migration
		fb=ceil(de);//simplest case: heaviside step function
	}else{//host migration
		fb=de;//simplest case with migration
	}

	return fb;
}
/****************************************
*	moviment of hosts:		*
*	2 neighbors exchange places	*
*****************************************/
void hostMoviment(int id1,int id2){
	int i,htmp;
	double bac_tmp;

	for(i=0; i<TYPES; ++i){
		bac_tmp=bac[id1][i];
		bac[id1][i]=bac[id2][i];
		bac[id2][i]=bac_tmp;
	}
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

	dt=adjustTimeStep(event->cprobE[event->usizeE-1]);
	dnumsteps=ceil(Dt_ref/dt);
	dt=Dt_ref/(double)dnumsteps;
	event->dtE=dt;

	return dnumsteps;
}
/****************************************************************
*	Host Layer Evolution (for a time interval=Dt_ref) for	*
*	the complete graph version				*
*****************************************************************/
int evolveHostCG(int dnumsteps,Event *event, TimeMeasures *meas){
	int i,idh,idlist,ide,idk,ne,nh,nb,nd,*listb,*listd;
	double nr;
	DynList empty_viz;
	
	nh=listh->usize;
	
	empty_viz.size=N-nh;
	empty_viz.usize=0;
	empty_viz.vec = (int *)calloc(empty_viz.size,sizeof(int));
	for(i=0; i<empty_viz.size; ++i){
		ide=listh->vec[i+nh];
		//if(host[ide]==0){
		empty_viz.vec[i]=listh->vec[i+nh];
		++empty_viz.usize;
	//	}
	}
                
	listb=(int *)calloc(dnumsteps,sizeof(int));
	listd=(int *)calloc(dnumsteps,sizeof(int));
	nb=0;
	nd=0;

	for(i=0; i<dnumsteps; ++i){
		nr=FRANDOM;
		if(nr<event->cprobE[event->usizeE-1]){
			event->whichE=selectEventCP(nr,event->cprobE,event->usizeE);
			idh=listh->vec[event->whichE%nh];
			if(host[idh]==1){//if chosen host is alive and the system has more than 1 host
				switch(event->whichE/nh){
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
							timeb[idk]=meas->Tnow+Dt_ref;//time of birth of host @idk
							if(timeb[idh]>=0.){
								meas->timegh+=(meas->Tnow+Dt_ref-timeb[idh]);//host @idh first reproduction time is the current time minos the time of its birth
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
	return dnumsteps;
}
/****************************************************************
*	Host Layer Evolution (for a time interval=Dt_ref) for	*
*	the lattice version					*
*****************************************************************/
void evolveHostSL(int dnumsteps,Event *event, TimeMeasures *meas){
	int i,idh,idlist,idlistk,ide,idk,ne,nh,nb,nd,nm,*listb,*listd,*listm;
	double nr1,nr2,fb,fm;
	double mh=spar->migh;
	DynList empty_viz;
	
	nh=listh->usize;
	empty_viz.usize=0;
	empty_viz.size=VIZ;
	empty_viz.vec = (int *)calloc(empty_viz.size,sizeof(int));
                
	listb=(int *)calloc(dnumsteps,sizeof(int));
	listd=(int *)calloc(dnumsteps,sizeof(int));
	nb=0;
	nd=0;
	if(mh>0.){
		listm=(int *)calloc(dnumsteps,sizeof(int));
	}
	nm=0;

	for(i=0; i<dnumsteps; ++i){
		nr1=FRANDOM;
		if(nr1<event->cprobE[event->usizeE-1]){
			event->whichE=selectEventCP(nr1,event->cprobE,event->usizeE);
			idh=listh->vec[event->whichE%nh];
			if(host[idh]==1){//if chosen host is alive and the system has more than 1 host
				switch(event->whichE/nh){
					case 0: 
						searchEmptyNeighbors(0,idh,host,neighbor,&empty_viz);
						ne=empty_viz.usize;//number of empty sites in the neighborhood

						fb=birthFunc(ne,VIZ,BFunc);/*if there is no host migration, returns heaviside(ne), ne/VIZ otherwise*/
						nr2=FRANDOM;
						fm=(1.-fb)*mh;
						if(nr2<fb){//birth
							ide=randNeighborID(idh,empty_viz.vec,0,ne,empty_viz.size);//randomly chooses the index of an empty site, stored on empty_viz.vec
							idk=empty_viz.vec[ide];//host id of the empty site, idk, is the ide-th element of empty_viz.vec
							listb[nb]=idk;//list of sites that receive offspring during the Dt_ref time interval
							++nb;
							hostBirth(idh,idk);
							//measures	
							#ifdef GENTIME
							timeb[idk]=meas->Tnow+Dt_ref;//time of birth of host @idk
							if(timeb[idh]>=0.){
								meas->timegh+=(meas->Tnow+Dt_ref-timeb[idh]);//host @idh first reproduction time is the current time minos the time of its birth
								timeb[idh]=-1.;//just the first reproduction of an individual counts
								++meas->ngh;
							}
							#endif
						}else if(nr2<fm+fb){//host migration
							listm[nm]=idh;
							++nm;	
						}
						break;
					case 1: 
						listd[nd]=idh;//list of sites that died during the Dt_ref time interval
						++nd;
						hostDeath(idh);
						#ifdef GENTIME
						timeb[idh]=-1.;
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
		exchange(inverselisth,listh->vec[nh],idh);
		exchange(listh->vec,nh,idlist);
		++listh->usize;
	}
	for(i=0; i<nd; ++i){//death list
		nh=listh->usize;
		idh=listd[i];
		idlist=inverselisth[idh];
		exchange(inverselisth,idh,listh->vec[nh-1]);
		exchange(listh->vec,idlist,nh-1);
		host[idh]=0;
		--listh->usize;
	}
	for(i=0; i<nm; ++i){//host migration list
		idh=listm[i];
		if(host[idh]==1){//host is still alive
			idk=neighbor[idh][(int)(FRANDOM*VIZ)];//random neighbor
			hostMoviment(idh,idk);//idh and idviz change places
			
			idlist=inverselisth[idh];
			idlistk=inverselisth[idk];
			exchange(inverselisth,idh,idk);
			exchange(listh->vec,idlist,idlistk);

		}
	}
	
	/**************/
	#ifdef NUMHEVENTSxT
	meas->numb=nb;
	meas->numd=nd;
	#endif

	if(spar->migh>0.){
		free(listm);
	}
	free(empty_viz.vec);
	free(listb);
	free(listd);
	return;
}
