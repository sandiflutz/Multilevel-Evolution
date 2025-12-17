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
void calcHostEvents(double *event){
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
		w[i]=calcAcumInvest(idh);
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
        double fk,fp,mean,*cprob=NULL;


        cprob=(double *)calloc(TYPES,sizeof(double));//cumulative probability for the bacteria types
	cprob[0]=bac[idp][0]/spar->micr[idp];
        for(i=1; i<TYPES; ++i){
                cprob[i]=cprob[i-1]+bac[idp][i]/spar->micr[idp];
        }

        ns=0;
	mean=10e03*Bacv;
	sample=poissonRandNum(mean);
        while(ns<sample){
		p=selectEventCP(cprob,TYPES);
                bac[idk][p]+=(double)Bacv/sample;
        }

	spar->micr[idk]=Bacv;

        free(cprob);
        return;
}
/******************************************************
 *      birth of a new host                           *
 *****************************************************/
void hostBirth(int idp, int idk){

        host[idk]=2;//newborn
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

        host[idh]=0;//newborn
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
        double passive_time;

        passive_time=exp(sumprob);

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
void dynamicsHost(Event *event){
        int id,idh,id_list,idk,nh;
	int idev=event->whichE;

	nh=listh->usize;//# of alive hosts

        id=idev%nh;/*If idev<nh, id=idev:birth event; otherwise nh<=idev<2nh, id=idev-nh:death event*/
        idh=listh->vec[id];

        if(idev<nh){//reproduction of host @idh
                #if (NETWORK==0)//well-mixed
                int ne=N-nh;
		if(ne>0){//there are empty sites in the system
                	id_list=(int)(FRANDOM*ne)+nh-1;//randomly select an index from the second part of the host list, where empty sites are stored 
              		idk=listh->vec[id_list];//position index of the empy site (that is going to receive the host offspring)
			hostBirth(idh,idk);
			//updating host's lists
			listAdd(listh,idk,id_list);
			listSimpleAdd(listnb,idk);
		}
                #else//square-lattice 
		searchLiveNeighbors(0,idh,VIZ,host,neighbors,alive_viz);
		int nav=alive_viz->usize;
                if(VIZ-nav>0){//there are empty sites in the neighborhood of host[idh]
                        id_list=(int)(FRANDOM*(VIZ-aviz))+nav-1;
                        idk=alive_viz.vec[id_list];
			hostBirth(idh,idk);
                	//updating host's lists
			listAdd(listh,idk,id_list);
			listSimpleAdd(listnb,idk);
                }
                #endif
        }else{//death of host @idh
		hostDeath(idh);
		listSub(listh,idh,id);/*exchange elements @idh, position @id (which is in the fisrt part of the list of 
					alive hosts: elements from 0 to @nh-1 on the vector that is part of the struct @listh), 
					with the (@nh-1)-ith element, then decrement @nh
					*/
        }

        return;
}
