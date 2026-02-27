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
	double cinv,kbac=spar->kbac;


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
/***************************************************************
*     Set microbial frequencies for the offspring of hosr @idp *
*     Bacteria types and their frequencies are randomly        *
*     selected using a normal distrution around the bacteria   *
*     frequencies on the parent                                *
****************************************************************/
void setMicrKidsNorm(int idp, int idk){
        int i,ns,p;
        double fk,fp,norm,*cprob=NULL;


        cprob=(double *)calloc(TYPES,sizeof(double));//cumulative probability for the bacteria types
	cprob[0]=bac[idp][0]/spar->micr[idp];
        for(i=1; i<TYPES; ++i){
                cprob[i]=cprob[i-1]+bac[idp][i]/spar->micr[idp];
        }

        norm=0.;
        ns=0;
        while(ns<BSAMPLES){
		p=selectEventCP(cprob,TYPES);
		fp=bac[idp][p]/spar->micr[idp];
                fk=truncGaussRandNum(fp,spar->sigma,0.,1.);
                bac[idk][p]+=fk;
                norm+=fk;
                ++ns;
        }

        for(i=0; i<TYPES; ++i){
                bac[idk][i]=bac[idk][i]*Bacv/norm;//normalizing bac so the sum of microbes in host idk is bv
        }

	spar->micr[idk]=Bacv;

        free(cprob);
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
        double mean,*cprob=NULL;


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
		p=selectEventCP(cprob,TYPES);
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
	#else//vertical transmission using a poisson dist.
	setMicrKidsPoiss(idp,idk);//set new host microbiome
	#endif


        return;
}
/******************************************************
*       death of a host                               *
*******************************************************/
void hostDeath(int idh){
        int i;

        host[idh]=0;
        for(i=0; i<TYPES; ++i){
                bac[idh][i]=0.;
        }
	spar->micr[idh]=0.;

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
/********************************************************
*                  host dynamics                        *
* Receives the id of the focus host, the type of event  *
* (reproduction=1, or death=0, the list with the labels *
* of existing hosts, lhost[host]=label, and the inverse *
* hosts list, ilhost[label]=host)                       *                 
*********************************************************/
void dynamicsHost(int type_event,int idh,DynList *lhost,int *ilhost){
        int idlist_h,idlist_k,idk,nh;
	
	idlist_h=ilhost[idh];//label of the focus host on the temporary list of alive hosts 
	nh=lhost->usize;
	switch(type_event){
		case 0://reproduction
			#if (NETWORK==0)//complete graph
			idlist_k=randNeighbor(idlist_h,lhost->vec,lhost->usize,N-1,N);/*choose a random neighbor from the host list. Variable being sent: 1-focus id,2-neighbors vector,
											 *3 and 4-randomly select neighbors between these 2 ids,
											 *5-neighbors vector size*/
                        idk=lhost->vec[idlist_k];			
			#else
			searchLiveNeighbors(0,idh,VIZ,host,neighbors,alive_viz);
			nav=alive_viz->usize;//# of alive neighbors
			idk=randNeighbor(idh,alive_viz->vec,nav,VIZ-1,VIZ);/*passing: 1-focus host id,2-list of alive neighboring hosts (including the focus host),
									      *3- position id of the first empty space = # of alive neigbors, 4-position id of the last empty space
									      *5-size of the list of neighbors = maximum number of alive neighbors
									      *the last element of @alive_viz, the (VIZ-1)th element, 
									      *is the focus host, so it cannot be include. */
			idlist_k=ilhost[idk];
			#endif
			exchange(ilhost,lhost->vec[nh],idk);
			exchange(lhost->vec,nh,idlist_k);
			++lhost->usize;
			hostBirth(idh,idk);
			break;
		case 1: //host death
			hostDeath(idh);
			exchange(ilhost,idh,lhost->vec[nh-1]);
			exchange(lhost->vec,idlist_h,nh-1);
			--lhost->usize;
			break;
	}


        return;
}
/*************************************************************
*    Host Layer Evolution (for a time interval=Dt_ref):      *
*    This version uses a dynamical value for the time step   *
*    of the host layer of dt<=Dt_ref                         *
*    (original paper version)                                *
**************************************************************/
int evolveHostDtH(Event *event,DynList *listh_tmp, int *inverselisth_tmp,TimeMeasures *meas){
	int i,idh,nh,dnumsteps;
	double dt,sumprobs;
	                
	dt=adjustTimeStep(event->cprobE[event->usizeE-1]);
	dnumsteps=ceil(Dt_ref/dt);
	dt=Dt_ref/(double)dnumsteps;
	sumprobs=dt*event->cprobE[event->usizeE-1];
	event->dtE=dt;
                
	nh=listh->usize;
	/*temporary host list (to keep track of the changes in the host list, that have to be updated after the host time substeps)*/
	listh_tmp->usize=listh->usize;
	for(i=0; i<N; ++i){
		listh_tmp->vec[i]=listh->vec[i];
		inverselisth_tmp[i]=inverselisth[i];
	}
	/******/

	for(i=0; i<dnumsteps; ++i){
		if(FRANDOM<sumprobs){
			event->whichE=selectEventCP(event->cprobE,event->usizeE);
			idh=listh->vec[event->whichE%nh];
			if((host[idh]==1)&&(listh_tmp->usize>1)){//if chosen host is alive and the system has more than 1 host
				switch(event->whichE/nh){
					case 0: dynamicsHost(0,idh,listh_tmp,inverselisth_tmp);
						break;
					case 1: dynamicsHost(1,idh,listh_tmp,inverselisth_tmp);
						break;
				}
				#ifdef NUMHEVENTSxT
				meas->numb+=1-event->whichE/nh;
				meas->numd+=event->whichE/nh;
				#endif
			}
		}
	}
	                
	/*updating hosts dynamic list*/
	for(i=0; i<N; ++i){
		listh->vec[i]=listh_tmp->vec[i];
		inverselisth[i]=inverselisth_tmp[i];
	}
	listh->usize=listh_tmp->usize;
	/**************/

	return dnumsteps;
}
/*************************************************************
*    Host Layer Evolution (for a time interval=Dt_ref):      *
*    This version uses a tau-leaping method                  *
**************************************************************/
void evolveHostTLP(Event *event,DynList *listh_tmp, int *inverselisth_tmp,TimeMeasures *meas){
	int i,k,idh,nh,evorder;                
	double pb;

	nh=listh->usize;
	/*temporary host list (to keep track of the changes in the host list, that have to be updated after the host time substeps)*/
	listh_tmp->usize=listh->usize;
	for(i=0; i<N; ++i){
		listh_tmp->vec[i]=listh->vec[i];
		inverselisth_tmp[i]=inverselisth[i];
	}
	/***/
		
	k=poissonRandKnuth((event->dtE*event->cprobE[event->usizeE-1]));
	for(i=0; i<k; ++i){
		event->whichE=selectEventCP(event->cprobE,event->usizeE);
		idh=listh->vec[event->whichE%nh];
		if((host[idh]==1)&&(listh_tmp->usize>1)){//if chosen host is alive and the system has more than 1 host
			switch(event->whichE/nh){
				case 0: dynamicsHost(0,idh,listh_tmp,inverselisth_tmp);
					break;
				case 1: dynamicsHost(1,idh,listh_tmp,inverselisth_tmp);
					break;
			}
			#ifdef NUMHEVENTSxT
			meas->numb+=1-event->whichE/nh;
			meas->numd+=event->whichE/nh;
			#endif
		}
	}
		
	//updating hosts dynamic list
	for(i=0; i<N; ++i){
		listh->vec[i]=listh_tmp->vec[i];
		inverselisth[i]=inverselisth_tmp[i];
	}
	listh->usize=listh_tmp->usize;
	nh=listh->usize;
	/**************/
	return;
}
