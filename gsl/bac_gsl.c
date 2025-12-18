/* bac_gsl.c */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include"randgen_ufrgs.h"
#include"globals.h"
#include"tools.h"
#include"bac_gsl.h"
/********************************************
*     GSL-compatible RHS function           *
*********************************************/
int bac_rhs(double t, const double y[], double f[], void *params){
	int j,jpl,jmi;
	double birth,death;
	SysParams *sp = (SysParams*) params;
	int idh=sp->idhost;
	double mu=sp->mu;//mutation rate
	double cost=sp->cost;//cost of helping for an ideal helper
	double beta=sp->beta;//birth rate of a neutral bacteria
	double delta=sp->delta;//bacteria death rate
	
	memset(f,0.,sizeof(double)*TYPES);

	for(j=0; j<TYPES;++j){
		/*births*/
		birth=(1.-mu)*(1.-cost*sp->s[j]*sp->inv[j])*beta*y[j];//division of type j
		jpl=j+1;
		jmi=j-1;
		if(j==0){
			jmi=jpl;
		}else if(j==TYPES-1){
			jpl=jmi;
		}
		birth+=0.5*mu*(1.-cost*sp->s[jpl]*sp->inv[jpl])*beta*y[jpl];//division of type j+1 -> mutation into j
		birth+=0.5*mu*(1.-cost*sp->s[jmi]*sp->inv[jmi])*beta*y[jmi];//division of type j-1 -> mutation into j
		/*death*/
		death=delta*sp->micr[idh]*y[j];
		/*bacteria function without the migration terms*/
		f[j]=birth-death;
	}

	return GSL_SUCCESS;
}
/********************************************
*     GSL-compatible jacobian               *
*********************************************/
int bac_jac(double t, const double y[], double *dfdy,double *dfdt, void *params){
	int b,bpl,bmi,idj1,idj2,ljac;
	double diag_par;
	SysParams *sp = (SysParams*) params;
	int idh=sp->idhost;//index of the host related to this jacobian
	double mu=sp->mu;//mutation rate
	double cost=sp->cost;//cost of helping for an ideal helper
	double beta=sp->beta;//birth rate of a neutral bacteria
	double delta=sp->delta;//bacteria death rate

	ljac=TYPES;//linear size of the jacobian matrix
	
	memset(dfdt,0.,sizeof(double)*ljac);	
	memset(dfdy,0.,sizeof(double)*ljac*ljac);

	/*the only nonzero elements of the jacobian of the bacteria function (without the migration terms) 
	 * are the diagonal and the sub and supra diagonal (because of the mutation terms type->type+1, type->type-1*/

	for(b=0; b<TYPES; ++b){
		/*diagonal terms:same type for f[] and y[] in the df/dy[] */
		diag_par=(1.-mu)*(1.-cost*sp->s[b]*sp->inv[b])*beta-delta*sp->micr[idh];
		
		idj1=b*TYPES+b;//diagonal index of a linearized jacobian matrix of size TYPESxTYPES
		dfdy[idj1]=diag_par-y[b];//diagonal jac term
					     
		/*sub and supra diagonals: f[b] being derived by types b+1 and b-1*/
		bpl=b+1;
		idj1=b*TYPES+bpl;
	
		bmi=b-1;
		idj2=b*TYPES+bmi;
		if(b==0){//if type=0, mutations happen only from type->type+1
			bmi=bpl;  
			idj2=idj1;
		}else if(b==TYPES-1){//if type=TYPES-1, mutations happen only from type->type-1
			bpl=bmi;
			idj1=idj2;
		}

		dfdy[idj1]=0.;
		dfdy[idj2]=0.;
	
		dfdy[idj1]+=0.5*mu*beta*(1.-cost*sp->s[bpl]*sp->inv[bpl]);//supradiagonal jac term
		
		dfdy[idj2]+=0.5*mu*beta*(1.-cost*sp->s[bmi]*sp->inv[bmi]);//subdiagonal jac term
			
	}


	return GSL_SUCCESS;
}
/********************************************
*     Buil a GSL  ODE system struct         *
*********************************************/
void bac_make_system(gsl_odeiv2_system *sys, SysParams *P){
	sys->function = bac_rhs;
	sys->jacobian = bac_jac;
	sys->dimension = TYPES;
	sys->params=P;

	return;
}
/**************************************************
 *          Bacteria Dynamics                     *
 **************************************************/
int bacDynamics(gsl_odeiv2_driver *driver,Event *event){
	int i,j,k,idh,idk,nh,status;
	double migr_in,migr_out,time_tmp,**bac_tmp=NULL;
	double dt=event->dtE;
        double mig=spar->mig;//migration rate


	nh=listh->usize;//# of alive hosts
	if(nh>1){
		//allocating space for a tempory bacteria abundance vector (so the abundances can be updated synchronously)
		bac_tmp=(double **)calloc(N,sizeof(double *));
		for(i=0; i<N; ++i){
			bac_tmp[i]=(double *)calloc(TYPES,sizeof(double));
		}
		/*************/
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];//index of the i-th alive host stored in a list
			spar->idhost=idh;
			#if (NETWORK!=0)//not the complete graph
			searchLiveNeighbors(0,idh,host,neighbor,alive_viz);//store the indexes of @idh alive neighbors in a list
			int nv=alive_viz->usize;//# of alive neighbors
			#endif
			spar->micr[idh]=0.;
			for(j=0; j<TYPES; ++j){
				/*migration terms*/
				migr_in=0.;
				migr_out=mig*bac[idh][j];	
				#if (NETWORK==0)//well-mixed
				for(k=0; k<nh; ++k){
					idk=listh->vec[k];
					migr_in+=mig*bac[idk][j];
				}
				migr_in=(migr_in-mig*bac[idh][j])/nh;
				#else
				for(k=0; k<nv; ++k){
					idk=alive_viz->vec[k];
					migr_in+=mig*bac[idk][j];
				}
				migr_in/=nv;
				#endif
				//evolution due just to the migrations terms
				bac_tmp[idh][j]=bac[idh][j]+(migr_in-migr_out)*dt;
				spar->micr[idh]+=bac_tmp[idh][j];
			}
			gsl_odeiv2_driver_reset(driver);
			time_tmp=event->timeE;
			status=gsl_odeiv2_driver_apply(driver, &time_tmp, time_tmp+dt, bac_tmp[idh]);//evolution due to birth and death for host @idh
			if (status != GSL_SUCCESS) {
				for(j=0; j<N; ++j){
					free(bac_tmp[j]);
				}
				free(bac_tmp);
				printf("Error in ODE integration: %s\n", gsl_strerror(status));
				return 1;
			}

	
		}
	
		/*synchronously updating hosts microbiome abundances and bacteria abundances*/
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			spar->micr[idh]=0.;
			for(j=0; j<TYPES; ++j){
				bac[idh][j]=bac_tmp[idh][j];
				spar->micr[idh]+=bac[idh][j];
			}
		}
	
		/*freeing temporary array for bacteria state*/
		for(i=0; i<N; ++i){
			free(bac_tmp[i]);
		}
		free(bac_tmp);
	}else{
		idh=listh->vec[0];
		spar->idhost=idh;
		gsl_odeiv2_driver_reset(driver);
		time_tmp=event->timeE;
		status=gsl_odeiv2_driver_apply(driver, &time_tmp, time_tmp+dt, bac[idh]);//evolution due to birth and death for host @idh
		if (status != GSL_SUCCESS) {
			printf("Error in ODE integration: %s\n", gsl_strerror(status));
			return 1;
			
		}
		spar->micr[idh]=0.;
		for(j=0; j<TYPES; ++j){
			spar->micr[idh]+=bac[idh][j];
		}
	}


	return 0;
}
