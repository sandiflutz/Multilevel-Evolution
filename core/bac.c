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
* 	using euler method when the cost of helping	*
* 	depends bacteria abundances			*
 ********************************************************/
void bacDynCostVec(int idh,double dt,double *btmp,double time,double *type_abund,int nv){
	int j,jpl,jmi;
	double birth,death,migr_in,migr_out,func;
	double mut=spar->mu;//mutation rate
        double birthr=Beta;//birth rate of a neutral bacteria
        double deathr=Delta;//bacteria death rate
        double mig=spar->mig;//migration rate
	
	for(j=0; j<TYPES;++j){
		/*births*/
		birth=(1.-mut)*(1.-costvec[idh][j])*birthr*bac[idh*TYPES+j];//division of type j
		
		jpl=j+1;
		jmi=j-1;
		if(jpl<TYPES){
			birth+=0.5*mut*(1.-costvec[idh][jpl])*birthr*bac[idh*TYPES+jpl];//division of type j+1 -> mutation into j
		}
		if(jmi>=0){
			birth+=0.5*mut*(1.-costvec[idh][jmi])*birthr*bac[idh*TYPES+jmi];//division of type j-1 -> mutation into j
		}
		/*death*/
		death=deathr*spar->micr[idh]*bac[idh*TYPES+j];
		/*migrations*/
		migr_out=mig*bac[idh*TYPES+j];//emmigration

		//immigration
		migr_in=0.;
		if(nv>1){//nv is the number of alive hosts in the neighborhood of idh, including idh
			migr_in=mig*(type_abund[j]-bac[idh*TYPES+j])/((double)nv-1.);//type_abund[j]=abundance of type j in the neighborhood=sum_viz bac[viz][j] (where viz includes idh);
		}
		func=birth-death-migr_out+migr_in;
		btmp[idh*TYPES+j]=bac[idh*TYPES+j]+func*dt;
	}
	
	return;
}
/********************************************************
* 	integration of the microbial equations		*
* 	using euler method				*
 ********************************************************/
void bacDyn(int idh,double dt,double *btmp,double time,double *type_abund,int nv){
	int j,jpl,jmi;
	double birth,death,migr_in,migr_out,func;
	double mut=spar->mu;//mutation rate
	double cost=spar->cost;//cost of helping for the most helpfull type
        double birthr=Beta;//birth rate of a neutral bacteria
        double deathr=Delta;//bacteria death rate
        double mig=spar->mig;//migration rate


	for(j=0; j<TYPES;++j){
		/*births*/
		birth=(1.-mut)*(1.-cost*spar->inv[j])*birthr*bac[idh*TYPES+j];//division of type j
		
		jpl=j+1;
		jmi=j-1;
		if(jpl<TYPES){
			birth+=0.5*mut*(1.-cost*spar->inv[jpl])*birthr*bac[idh*TYPES+jpl];//division of type j+1 -> mutation into j
		}
		if(jmi>=0){
			birth+=0.5*mut*(1.-cost*spar->inv[jmi])*birthr*bac[idh*TYPES+jmi];//division of type j-1 -> mutation into j
		}
		/*death*/
		death=deathr*spar->micr[idh]*bac[idh*TYPES+j];

		//horizontal transmission
		migr_in=0.;
		migr_out=0.;
		if(nv>1){//nv is the number of alive hosts in the neighborhood of idh, including idh
			migr_out=mig*bac[idh*TYPES+j];//emmigration
			migr_in=mig*(type_abund[j]-bac[idh*TYPES+j])/((double)nv-1.);//type_abund[j]=abundance of type j in the neighborhood=sum_viz bac[viz][j] (where viz includes idh);
		}
		func=birth-death-migr_out+migr_in;
		btmp[idh*TYPES+j]=bac[idh*TYPES+j]+func*dt;
        }

	return;
}
/****************************************************************
 * Microbial layer evolution (for a time interval=dt)	*
 ****************************************************************/
void evoBac(double dt,double time){
	int i,j,idh,nh;
	double *bac_tmp,*type_abund;

	nh=listh->usize;
	type_abund=(double *)calloc(TYPES,sizeof(double));
	memset(type_abund,0.,sizeof(double)*TYPES);

	#if (NETWORK==0)//complete graph
	for(j=0; j<TYPES; ++j){
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			type_abund[j]+=bac[idh*TYPES+j];
		}
	}
	#else
	int nv,k,idviz;
	#endif

	
	bac_tmp=(double *)calloc(N*TYPES,sizeof(double));
	
	#if (Tneg>0)//with negative type the cost depends on negative types abundances
        for(i=0; i<nh; ++i){
                idh=listh->vec[i];
		updateCosts(idh,spar->micr[idh],spar->inv,spar->cost,costvec,bac);
	}
        for(i=0; i<nh; ++i){
                idh=listh->vec[i];

		#if (NETWORK==0)//complete graph
		bacDynCostVec(idh,dt,bac_tmp,time,type_abund,nh);
                #elif (NETWORK==1)//square lattice
		nv=(1.-rho_e[idh])*(double)VIZ+1;
		
		for(j=0; j<TYPES; ++j){
			type_abund[j]=bac[idh*TYPES+j];
		}

		if(nv>1){
			for(k=0; k<VIZ; ++k){
				idviz=neighbor[idh][k];
				if(host[idviz]==1){
					for(j=0; j<TYPES; ++j){
						type_abund[j]+=bac[idviz*TYPES+j];
					}
				}
			}
		}
		bacDynCostVec(idh,dt,bac_tmp,time,type_abund,nv);
                #else//small-world
		nv=(1.-rho_e[idh])*(double)con[idh]+1;
		
		for(j=0; j<TYPES; ++j){
			type_abund[j]=bac[idh*TYPES+j];
		}

		if(nv>1){
			for(k=0; k<con[idh]; ++k){
				idviz=neighbor[idh][k];
				if(host[idviz]==1){
					for(j=0; j<TYPES; ++j){
						type_abund[j]+=bac[idviz*TYPES+j];
					}
				}
			}
		}
		bacDynCostVec(idh,dt,bac_tmp,time,type_abund,nv);
                #endif
	}
	#else//no negative type: cost rate is gamma*inv[type]
        for(i=0; i<nh; ++i){
                idh=listh->vec[i];

		#if (NETWORK==0)//complete graph
		bacDyn(idh,dt,bac_tmp,time,type_abund,nh);
                #elif(NETWORK==1)//square lattice
		nv=(1.-rho_e[idh])*(double)VIZ+1;
		
		for(j=0; j<TYPES; ++j){
			type_abund[j]=bac[idh*TYPES+j];
		}

		if(nv>1){
			for(k=0; k<VIZ; ++k){
				idviz=neighbor[idh][k];
				if(host[idviz]==1){
					for(j=0; j<TYPES; ++j){
						type_abund[j]+=bac[idviz*TYPES+j];
					}
				}
			}
		}
		bacDyn(idh,dt,bac_tmp,time,type_abund,nv);
                #else//small-world
		nv=(1.-rho_e[idh])*(double)con[idh]+1;
		
		for(j=0; j<TYPES; ++j){
			type_abund[j]=bac[idh*TYPES+j];
		}

		if(nv>1){
			for(k=0; k<con[idh]; ++k){
				idviz=neighbor[idh][k];
				if(host[idviz]==1){
					for(j=0; j<TYPES; ++j){
						type_abund[j]+=bac[idviz*TYPES+j];
					}
				}
			}
		}
		bacDyn(idh,dt,bac_tmp,time,type_abund,nv);
                #endif
	}
	#endif



	//updating microbial states
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		spar->micr[idh]=0.;
		for(j=0;j<TYPES; ++j){
			bac[idh*TYPES+j]=bac_tmp[idh*TYPES+j];
			spar->micr[idh]+=bac[idh*TYPES+j];
		}
	}

	free(bac_tmp);
	free(type_abund);
	return;
}
