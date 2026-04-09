/* bac.c */
#include<math.h>
#include"randgen_ufrgs.h"
#include"globals.h"
#include"tools.h"
#include"bac.h"

/****************************************************************
 *     cost reduction function caused by negative types         *
 *     of bacteria                                              *
 ***************************************************************/
double expCostReduc(double x, double a,double b,double eps){
	double func;

	if(1.-x<=eps){
		func=a;
	}else{
		func=a*exp(-b*x/(1.-x));
	}
	
	return func;
}
/************************************************************************
 *	update costs multiplicative factor acording to the amount	* 
 *	of negative microbial types					*
 ************************************************************************/
void updateCosts(int idh,double mtot,double c0,double *w,double **cvec,double **bstate){
	int j;
	double mneg,eps=1e-6;
		
	mneg=0.;
		
	for(j=0; j<Tneg; ++j){
		mneg+=bstate[idh][j];
	}
	mneg/=mtot;

	for(j=0; j<Tneg; ++j){
		cvec[idh][j]=expCostReduc(mneg,CRnn0,CRnn1,eps);
		cvec[idh][j]*=c0*w[j];
	}
	for(j=Tneg; j<=Tpos+Tneg; ++j){
		cvec[idh][j]=expCostReduc(mneg,CRnp0,CRnp1,eps);
		cvec[idh][j]*=c0*w[j];
	}

	return;
}
/********************************************************
* 	integration of the microbial equations		*
* 	using euler method				*
 ********************************************************/
void bacDyn(int idh,double dt,DynList *liveviz,double **btmp,double time){
	int j,jpl,jmi,k,idviz,nv;
	double birth,death,migr_in,migr_out,func;
	
	double mut=spar->mu;//mutation rate
        double birthr=Beta;//birth rate of a neutral bacteria
        double deathr=Delta;//bacteria death rate
        double mig=spar->mig;//migration rate


	nv=liveviz->usize;
	for(j=0; j<TYPES;++j){
		/*births*/
		birth=(1.-mut)*(1.-costvec[idh][j])*birthr*bac[idh][j];//division of type j
		
		jpl=j+1;
		jmi=j-1;
		if(jpl<TYPES){
			birth+=0.5*mut*(1.-costvec[idh][jpl])*birthr*bac[idh][jpl];//division of type j+1 -> mutation into j
		}
		if(jmi>=0){
			birth+=0.5*mut*(1.-costvec[idh][jmi])*birthr*bac[idh][jmi];//division of type j-1 -> mutation into j
		}
		/*death*/
		death=deathr*spar->micr[idh]*bac[idh][j];
		/*migrations*/
		migr_out=mig*bac[idh][j];//emmigration

		//immigration
		migr_in=0.;
		for(k=0; k<nv; ++k){//the last element of the alive neighbors list @alive_viz.vec[nv-1] is the focus host id
			idviz=liveviz->vec[k];
			migr_in+=mig*bac[idviz][j];
		}
		if(nv>1){
			migr_in=(migr_in-mig*bac[idh][j])/((double)nv-1.);
		}
		func=birth-death-migr_out+migr_in;
		btmp[idh][j]=bac[idh][j]+func*dt;
        }

	return;
}
/****************************************************************
 * Microbial layer evolution (for a time interval=dt)	*
 ****************************************************************/
void evoBac(double dt,double time){
	int i,j,idh,nh;
	double **bac_tmp;
	
	bac_tmp=(double **)calloc(N,sizeof(double *));
	for(i=0; i<N; ++i){
		bac_tmp[i]=(double *)calloc(TYPES, sizeof(double));
	}
        
	nh=listh->usize;
	#if (Tneg>0)
        for(i=0; i<nh; ++i){
                idh=listh->vec[i];
		updateCosts(idh,spar->micr[idh],spar->inv,spar->cost,costvec,bac);
	}
	#endif

        for(i=0; i<nh; ++i){
                idh=listh->vec[i];

		#if (NETWORK==0)//complete graph
		bacDyn(idh,dt,listh,bac_tmp,time);
                #else//square lattice
		searchLiveNeighbors(1,idh,host,neighbor,alive_viz);
		bacDyn(idh,dt,alive_viz,bac_tmp,time);
                #endif
	}
	



	//updating microbial states
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		spar->micr[idh]=0.;
		for(j=0;j<TYPES; ++j){
			bac[idh][j]=bac_tmp[idh][j];
			spar->micr[idh]+=bac[idh][j];
		}
	}

	for(i=0; i<N; ++i){
		free(bac_tmp[i]);
	}
	free(bac_tmp);
	return;
}
