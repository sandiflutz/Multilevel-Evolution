/* evo.c */

#include<math.h>
#include"globals.h"
#include"randgen_ufrgs.h"
#include"tools.h"
#include"evo.h"

/****************************************************************
*     calculates the accumulated investment in host @index      *
*****************************************************************/
double calcAcumInvest(int index,SysParams *spar){
        int j,nh;
	double cinv,kbac=spar->kbac;


        cinv=0.;
       	for(j=0; j<TYPES; ++j){
		cinv+=bac[index][j]*inv[j]/kbac;
        }

        return cinv;
}
/****************************************************
* Calculates the host event rates: birth and death  *
*****************************************************/
void calcHostEvents(double *event, SysParams *spar){
        int i,idh,idmax,nh;
	int gh=spar->gh;
	int kh=spar->kh;
	double beta=spar->beta;
	double sb=spar->sb;
	double sd=spar->sd;
	double *w=NULL;

	nh=listh->usize;

	w=(double *)calloc(nh,sizeof(double));

        idmax=2*nh;
        memset(event,0.,sizeof(double)*idmax);

	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		w[i]=calcAcumInvest(idh,spar);
	}

        //birth events
        if(nh<N){//there are empty site
                for(i=0; i<nh; ++i){
                        event[i]=(double)beta*(1.+sb*w[i])/gh;
                }
        }

        //death events
        if(nh>0){//there are alive hosts
                for(i=nh; i<idmax; ++i){
                        event[i]=beta*(1.-sd*w[i-nh])*nh/((double)kh*gh);
                }
        }

	free(w);
        return;
}
/***************************************************************
*     Set microbial frequencies of the offspring of a host     *
*     located. The parent and the kid hosts are located at     *
*     sites idp and idk, respectively                          *
****************************************************************/
void setMicrKids(int idp, int idk,SysParams *spar){
        int i,ns,p;
        double fk,fp,norm,*cprob=NULL;


        cprob=(double *)calloc(TYPES,sizeof(double));//cumulative probability for the bacteria types
	cprob[0]=bac[idp][0]/micr[idp];
        for(i=1; i<TYPES; ++i){
                cprob[i]=cprob[i-1]+bac[idp][i]/micr[idp];
        }

        norm=0.;
        ns=0;
        while(ns<BSAMPLES){
		p=selectEventCP(cprob,TYPES);
		fp=bac[idp][p]/micr[idp];
                fk=gaussRandNum(fp,spar->sigma,0.,1.);
                bac[idk][p]+=fk;
                norm+=fk;
                ++ns;
        }

        for(i=0; i<TYPES; ++i){
                bac[idk][i]=bac[idk][i]*Bacv/norm;//normalizing bac so the sum of microbes in host idk is bv
        }

	micr[idk]=Bacv;

        free(cprob);
        return;
}
/******************************************************
 *      birth of a new host                           *
 *****************************************************/
void hostBirth(int idp, int idk,SysParams *spar){

        host[idk]=2;//newborn
	setMicrKids(idp,idk,spar);//set new host microbiome


        return;
}
/******************************************************
*       death of a host                               *
*******************************************************/
void hostDeath(int idh){
        int i;

        host[idh]=0;//newborn
        for(i=0; i<TYPES; ++i){
                bac[idh][i]=0.;
        }
	micr[idh]=0.;

        return;
}
/************************************************
*         returns gillespie's time increment    *
************************************************/
double gillespieTime(double sumprob){
        double passive_time;

        passive_time=exp(sumprob);

        return passive_time;
}
/***********************************************
*         Adjust host time step                *
************************************************/
double adjustTimeStep(double *event,int size){
        int i,idm;
        double dt,dt_tmp,maxe1,maxe2,p2;
	DynListF *p;

	p=malloc(sizeof(DynListF));
	p->vecf=(double *)calloc(size,sizeof(double));
	p->sizef=size;
	p->usizef=size;

	for(i=0; i<size; ++i){
		p->vecf[i]=event[i];
	}

        dt=Dt_ref;

        idm=findMaxElement(p->usizef,&maxe1,p->vecf);
	listSimpleSubF(p,idm);
        
	idm=findMaxElement(p->usizef,&maxe2,p->vecf);
	if(maxe2==0.){
		if(maxe1>0.){
			dt_tmp=0.1/maxe1;
			if(dt>dt_tmp)dt=dt_tmp;
		}
	}else{
		p2=maxe1*maxe2;
		dt_tmp=0.1/sqrt(p2);
		if(dt>dt_tmp)dt=dt_tmp;
	}

	free(p->vecf);
	free(p);
        return dt;
}
/**************************************************
*      Adjust host time step: simple version      *
*      Chooses the maximum vector element @maxe   *
*      and considers that the maximum probability *
*      of 2 host events is ~maxe*maxe             * 
**************************************************/
double adjustTimeStepSimple(double *event,int size){
        int i,idm;
        double dt,dt_tmp,maxe;

        dt=Dt_ref;

        idm=findMaxElement(size,&maxe,event);
	if(maxe>0.){
        	dt_tmp=0.1/maxe;
        	if(dt>dt_tmp)dt=dt_tmp;
	}

        return dt;
}
/*******************************************************
*                  host dynamics                       *
********************************************************/
void dynamicsHost(Event *event,SysParams *spar){
        int id,idh,id_list,idk,ne,nav,nh;
	int idev=event->whichE;

	nh=listh->usize;//# of alive hosts

        id=idev%nh;/*If idev<nh, id=idev:birth event; otherwise nh<=idev<2nh, id=idev-nh:death event*/
        idh=listh->vec[id];

        if(idev<nh){//birth of host @list_host[k]
                #if (NETWORK==0)//well-mixed
                ne=N-nh;
		if(ne>0){//there are empty sites in the system
                	id_list=(int)(FRANDOM*ne)+nh-1;//randomly select an index from the second part of the host list, where empty sites are stored 
              		idk=listh->vec[id_list];//position index of the empy site (that is going to receive the host offspring)
			hostBirth(idh,idk,spar);
			//updating host's lists
			listAdd(listh,idk,id_list);
			listSimpleAdd(listnb,idk);
		}
                #else 
		searchLiveNeighbors(0,idh,VIZ,host,neighbors,alive_viz);
		nav=alive_viz->usize;
                if(VIZ-nav>0){//there are empty sites in the neighborhood of host[idh]
                        id_list=(int)(FRANDOM*(VIZ-aviz))+nav-1;
                        idk=alive_viz.vec[id_list];
			hostBirth(idh,idk,spar);
                	//updating host's lists
			listAdd(listh,idk,id_list);
			listSimpleAdd(listnb,idk);
                }
                #endif
        }else{//death of host @list_host[k-num_host]
		hostDeath(idh);
		listSub(listh,idh,id);/*exchange elements @idh, position @id (which is in the fisrt part of the list of 
					alive hosts: elements from 0 to @nh-1 on the vector that is part of the struct @listh), 
					with the (@nh-1)-ith element, then decrement @nh
					*/
        }

        return;
}
