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

	w=(double *)calloc(nh,sizeof(double));

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
/*******************************************************
 *  	calculate the number of times steps            *
 *******************************************************/
void calcNumSteps(int *dnt_h,int *dnt_b,int *dnt,double dt_h,double dt_b){
	
	*dnt_h=(int)(dt_b/dt_h);

	if(*dnt_h==0){
		*dnt_h = 1;
		*dnt_b = round(dt_h/dt_b);
		*dnt = *dnt_b;

	}else{
		*dnt_b=1;
		*dnt_h=round(dt_b/dt_h);
		*dnt = *dnt_h;
	}

	return;
}
/*******************************************************
*                  host dynamics                       *
********************************************************/
void dynamicsHost(int type_event,int idh,DynList *lhost,int *ilhost){
        int idlist_h,idlist_k,idk,nh;
	
	idlist_h=ilhost[idh];
	nh=lhost->usize;
	switch(type_event){
		case 0:
			#if (NETWORK==0)//complete graph
			idlist_k=randNeighbor(idlist_h,lhost->vec,lhost->usize,N-1,N);/*sending: 1-focus id,2-neighbors vector,
											 *3 and 4-randomly select neighbors between these 2 ids,
											 *5-neighbors vector size*/
                        idk=lhost->vec[idlist_k];			
			#else
			searchLiveNeighbors(0,idh,VIZ,host,neighbors,alive_viz);
			nav=alive_viz->usize;//# of alive neighbors
			idk=randNeighbor(idh,alive_viz->vec,nav,VIZ-1,VIZ);
			idlist_k=ilhost[idk];
			#endif
			exchange(ilhost,lhost->vec[nh],idk);
			exchange(lhost->vec,nh,idlist_k);
			++lhost->usize;
			hostBirth(idh,idk);
			break;
		case 1: hostDeath(idh);
			exchange(ilhost,idh,lhost->vec[nh-1]);
			exchange(lhost->vec,idlist_h,nh-1);
			--lhost->usize;
			break;
	}


        return;
}
