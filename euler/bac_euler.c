/* bac_euler.c */
#include"randgen_ufrgs.h"
#include"globals.h"
#include"tools.h"
#include"bac_euler.h"

void bac_euler(double dt,SysParams *sp){
	int i,j,jpl,jmi,idh1,k,idh2,nh;
	double birth,death,migr_in,migr_out,func,**bac_tmp;
	double mut=sp->mu;//mutation rate
        double cost=sp->cost;//cost of helping for an ideal helper
        double birthr=sp->beta;//birth rate of a neutral bacteria
        double deathr=sp->delta;//bacteria death rate
        double mig=sp->mig;//migration rate

	bac_tmp=(double **)calloc(N,sizeof(double *));
	for(i=0; i<N; ++i){
		bac_tmp[i]=(double *)calloc(TYPES, sizeof(double));
	}

        nh=listh->usize;
	
        for(i=0; i<nh; ++i){
                idh1=listh->vec[i];
                
		#if (NETWORK!=0)
                searchLiveNeighbors(0,idh1,host,neighbor,alive_viz);
		int nv=alive_viz->usize;
                #endif

                for(j=0; j<TYPES;++j){
                        /*births*/
                        birth=(1.-mut)*(1.-cost*sp->s[j]*sp->inv[j])*birthr*bac[idh1][j];//division of type j
                        jpl=j+1;
                        jmi=j-1;
                        if(j==0){
                                jmi=jpl;
                        }else if(j==TYPES-1){
                                jpl=jmi;
                        }
                        birth+=0.5*mut*(1.-cost*sp->s[jpl]*sp->inv[jpl])*birthr*bac[idh1][jpl];//division of type j+1 -> mutation into j
                        birth+=0.5*mut*(1.-cost*sp->s[jmi]*sp->inv[jmi])*birthr*bac[idh1][jmi];//division of type j-1 -> mutation into j
                        /*death*/
                        death=deathr*sp->micr[idh1]*bac[idh1][j];
                        /*migrations*/
                        migr_out=mig*bac[idh1][j];//emmigration

                        #if (NETWORK==0)//complete graph
			migr_in=0.;
                        for(k=0; k<nh;++k){
                                idh2=listh->vec[k];
                                migr_in+=mig*bac[idh2][j];
                        }
                        migr_in=(migr_in-mig*bac[idh1][j])/(nh-1.);
                        #else
			migr_in=0.;
			for(k=0; k<nv-1; ++k){//the last element of the alive neighbors list @alive_viz.vec[nv-1] is the focus host id
                                idh2=alive_viz.vec[k][idh1];
                                migr_in+=mig*bac[idh2][j];
                        }
                        migr_in/=(double)nv;
                        #endif

                        func=birth-death-migr_out+migr_in;

			bac_tmp[idh1][j]=bac[idh1][j]+func*dt;
                }
		
        }

	for(i=0; i<nh; ++i){
		idh1=listh->vec[i];
		sp->micr[idh1]=0.;
		for(j=0;j<TYPES; ++j){
			bac[idh1][j]=bac_tmp[idh1][j];
			sp->micr[idh1]+=bac[idh1][j];
		}
	}

	for(i=0; i<N; ++i){
		free(bac_tmp[i]);
	}
	free(bac_tmp);
	return;
}
