/* meaures.c */
#include<time.h>
#include<math.h>
#include"globals.h"
#include"tools.h"
#include"measures.h"
#include"init.h"
#include"evo.h"

/********************************************************************************************************
 * Organized in 3 blocks: 	1)routines related to for general memory allocating and files. 		*
 * 				2)routines to calculate dynamical properties of the system		*
 * 				3)routines to store data on files					*
 * 				4)general measuring routines that which measuring routine to call	*	
 ********************************************************************************************************/


/********************************************************************************************************
* 					Routines for: 							*				
* 					-allocating memory 						*
* 					-liberating memory						*
* 					-opening files							*
* 					-closing files							*
*********************************************************************************************************/
/***********************************************************
*     allocate memory for global arrays and structs        *
*     related to time measuraments                         *
************************************************************/
void allocateMemTM(void){
	int i;

        
	#ifdef DIFBACOMPxT
	offcomp = malloc(sizeof(DynVec));
        if (!offcomp) { perror("malloc"); exit(1);}
        offcomp->vecf=(double *)calloc(spar->kh,sizeof(double));//max. number of neighbor + focus host
	memset(offcomp->vecf,0.,sizeof(double)*spar->kh);
        offcomp->sizef=spar->kh;
        offcomp->usizef=0;
	#endif
	#ifdef BESTCLUSTER_TIMES
	bestwlisth = malloc(sizeof(DynVec));
	if (!bestwlisth) { perror("malloc"); exit(1);}
	bestwlisth->vec=(int *)calloc(spar->kh,sizeof(int));
	bestwlisth->size=spar->kh;
	bestwlisth->usize=0;
	maxclw=malloc(sizeof(ClusterFullID));
	if (!maxclw){ perror("malloc"); exit(1);}
	#endif

	/**generic file struct to use for storing data during the system evolution**/
	gfile = malloc(sizeof(GenFile));
        if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=300;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");

	/**file name components**/

        char *ngeral=(char *)calloc(200,sizeof(char));
        char *ntneg = (char *)calloc(50,sizeof(char));
	int npar=2,mul,nf;
        char nparam[npar][10];
        double expo,param[npar];
        param[0]=spar->mu;
        param[1]=spar->mig;

	for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
			expo=floor(log10(param[i]));
        		mul=ceil(param[i]/pow(10.,expo));
			nf=0;
			while(ceil(param[i]/pow(10.,expo))!=(int)(param[i]/pow(10.,expo))){
				--expo;
				mul=ceil(param[i]/pow(10.,expo));
				++nf;
			}
			sprintf(nparam[i],"%de%d",mul,(int)expo);
                }
        }
	#if (NETWORK==0)
		#if (CI!=1)
        	sprintf(ngeral,"CG_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_mh%0.1f_rmh%d_CI%d",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->mh,spar->rmigh,CI);
		#else
        	sprintf(ngeral,"CG_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_mh%0.1f_rmh%d",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->mh,spar->rmigh);
		#endif
	#elif (NETWORK==1)
		#if (CI!=1)
        	sprintf(ngeral,"SL_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_mh%0.1f_rmh%d_plr%0.2f_MT%d_CI%d",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->mh,spar->rmigh,spar->plr,MIGRATION_TYPE,CI);
		#else
        	sprintf(ngeral,"SL_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_mh%0.1f_rmh%d_plr%0.2f_MT%d",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->mh,spar->rmigh,spar->plr,MIGRATION_TYPE);
		#endif
	#else //smallworld
		#if (CI!=1)
        	sprintf(ngeral,"SW_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_mh%0.1f_rmh%d_psw%0.1f_CI%d",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->mh,spar->rmigh,Psw,CI);
		#else
        	sprintf(ngeral,"SW_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_mh%0.1f_rmh%d_psw%0.1f",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->mh,spar->rmigh,Psw);
		#endif
	#endif
                
	#if (Tneg>0)
	sprintf(ntneg,"_Tpos%d_Tneg%d_CRnnA%d_CRnnB%0.1f_CRnpA%d_CRnpB%0.1f",Tpos,Tneg,(int)CRnn0,CRnn1,(int)CRnp0,CRnp1);
	#else
	sprintf(ntneg,"");
	#endif
        sprintf(gfile->fname,"%s%s",ngeral,ntneg);
        free(ngeral);
        free(ntneg);

	return;
}
/************************************************
*       Free Allocated Memory for time		*
*       measurements				*
*************************************************/
void freeMemTM(void){
	
	if(gfile->file){
		fclose(gfile->file);
		gfile->file = NULL;
	}
	if(gfile->fname){
		free(gfile->fname);
		gfile->fname = NULL;
	}
	if(gfile->fdatapath){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile){
		free(gfile);
		gfile=NULL;
	}
	#ifdef DIFBACOMPxT
	if(offcomp->vecf){
		free(offcomp->vecf);
		offcomp->vecf=NULL;
	}
        if(offcomp){
		free(offcomp);
		offcomp=NULL;
	}
	#endif
	#ifdef BESTCLUSTER_TIMES
	if(bestwlisth->vec){
		free(bestwlisth->vec);
		bestwlisth->vec=NULL;
	}
	if(bestwlisth){
		free(bestwlisth);
		bestwlisth=NULL;
	}
	if(maxclw){
		free(maxclw);
		maxclw=NULL;
	}
	#endif

	return;
}
/******************************************
*	Open Global Files                 *
*******************************************/
void openFiles(void){

#ifdef MULTIPLE_COSTS_WxT
        int ok=0,namelen,dnl;
        unsigned long id;

        id = (unsigned long)time(NULL);

        dnl=200;
        namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;

        char *name=(char *)calloc(namelen,sizeof(char));
        
	while(ok==0){
                sprintf(name,"%saverInvXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
        sprintf(name,"%saverInvXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
        
	free(name);
#endif
#ifdef TMEAS
        int ok=0,namelen,dnl;
        unsigned long id;

        id = (unsigned long)time(NULL);

        dnl=200;
        namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;

        char *name=(char *)calloc(namelen,sizeof(char));

#ifdef AVERINVxT
        while(ok==0){
                sprintf(name,"%saverInvXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
        sprintf(name,"%saverInvXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef NUMHEVENTSxT
        while(ok==0){
                sprintf(name,"%snumheventsXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
        sprintf(name,"%snumheventsXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef CORRxT
        while(ok==0){
                sprintf(name,"%scorrXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
	sprintf(name,"%scorrXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef CLUSTERSxT
        while(ok==0){
                sprintf(name,"%sclustersXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
	sprintf(name,"%sclustersXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef CLUSTERS_DISTxT
        while(ok==0){
                sprintf(name,"%sclDistXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
	sprintf(name,"%sclDistXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef BESTCLUSTER_TIMES
        while(ok==0){
                sprintf(name,"%sclExtTimes_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
	sprintf(name,"%sclExtTimes_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef RVNxTxCORRBAC
        while(ok==0){
                sprintf(name,"%srvnXtXw_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
	sprintf(name,"%srvnXtXw_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef DIFBACOMPxT
        while(ok==0){
		#if (OFFCOMP==0)
                sprintf(name,"%soffdifXt_%s_Nb%d_%ld.dat",gfile->fdatapath,gfile->fname,SAMPLE,id);
		#else
                sprintf(name,"%soffainvXt_%s_Nb%d_%ld.dat",gfile->fdatapath,gfile->fname,SAMPLE,id);
		#endif
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
	#if (OFFCOMP==0)
	sprintf(name,"%savfinvpkdifXt_%s_Nb%d_%ld.dat",gfile->fdatapath,gfile->fname,SAMPLE,id);
	#else
	sprintf(name,"%soffainvXt_%s_Nb%d_%ld.dat",gfile->fdatapath,gfile->fname,SAMPLE,id);
	#endif
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
#endif

        free(name);
#endif
        return;
}
/**********************************************
*   		Close Global Files            *                    
***********************************************/
void closeFiles(void){
	
	if(gfile->file!=NULL){
		fclose(gfile->file);
		gfile->file=NULL;
	}

	return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************************
* 		Routines for calculating dynamical properties of the system			*
*************************************************************************************************/

/****************************************************************
 *	calculate investment density per Host:			*
 * 	densInvH[host]=sum_type(bac[host][type]*inv[type])    calcAverInv 	*
 ****************************************************************/
void calcInvDens(double *densInvH){
	int i,idh,j,nh;
	
	nh=listh->usize;
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		densInvH[i]=0.;
		for(j=0; j<TYPES; ++j){
			densInvH[i]+=bac[idh*TYPES+j]*spar->inv[j];
		}
	}

	return;
}
/************************************************************************
 * 	calculate investment frequency per Host:    			*
 * 	freqInvH[host]=sum_type(bac[host][type]*inv[type])/micr[host] 	*
 ************************************************************************/
void calcInvFreq(double *freqInvH){
	int i,idh,j,nh;
	
	nh=listh->usize;
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		freqInvH[i]=0.;
		for(j=0; j<TYPES; ++j){
			freqInvH[i]+=bac[idh*TYPES+j]*spar->inv[j]/spar->micr[idh];
		}
	}

	return;
}
/****************************************************************
 * 	calculate investment distribution among hosts		*
 ****************************************************************/
void calcInvDist(double binsize,double *hist_inv,double *freqInvH){
	int i,nh,id;

	nh=listh->usize;

	calcInvFreq(freqInvH);
	for(i=0; i<nh; ++i){
		id=(int)(freqInvH[i]/binsize);
		hist_inv[id]+=1.;
	}


	return;
}
/****************************************************************************************
* 	calculate current average investment in the system:				*
* 	averInv=sum_host(sum_type(bac[host][type]))/total amount of bac. in the system	*
****************************************************************************************/
double calcAverInv(void){
	int i,idh,nh;
	double tot_micr,averinv;
	double *densInvH=(double *)calloc(listh->usize,sizeof(double));

	nh=listh->usize;
	calcInvDens(densInvH);
	
	averinv=0.;
	tot_micr=0.;
        for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		tot_micr+=spar->micr[idh];
		averinv+=densInvH[i];
        }
	averinv/=tot_micr;//averaging over the wholcalcAverInve microbe population

	free(densInvH);	
	return averinv;
}
/****************************************************************
*	Calculate the average investment in each cluster, 	*
*	the label of the cluster with the highest 		*
*	investment and its size. Return the system average	*
*	investment.						*
*****************************************************************/
double findHighestInvCluster(int ncl,int *labels,double *wcl,double *clwstats,ClusterFullID *clid){
	int i,idh,nh;
	double w,averw,averwcl,averwcl2,ntot,*ntotcl;

	ntotcl=(double *)calloc(ncl,sizeof(double));

	nh=listh->usize;

	/*initializing the average investment of each cluster i, wcl[i], and the related total microbial abundance,ntotcl[i]*/
	for(i=0; i<ncl; ++i){
		wcl[i]=0.;
		ntotcl[i]=0.;
	}

	/*calculating wcl[i] and the average investment in the system*/
	averw=0.;
	ntot=0.;
	for(i=0; i<nh; ++i){		
		idh=listh->vec[i];
		w=calcAcumInvest(idh);
		wcl[labels[i]]+=w;
		averw+=w;
		ntotcl[labels[i]]+=spar->micr[idh];
		ntot+=spar->micr[idh];
	}
	averw/=ntot;
	clid->stateCL=0.;
	clid->whichLBF=0;
	averwcl=0.;
	averwcl2=0.;
	for(i=0; i<ncl; ++i){
		wcl[i]/=ntotcl[i];
		if(wcl[i]>clid->stateCL){
			clid->stateCL=wcl[i];
			clid->whichLBF=i;
		}
		averwcl+=wcl[i];
		averwcl2+=wcl[i]*wcl[i];
	}

	if(ncl>1){
		averwcl=(averwcl-wcl[clid->whichLBF])/((double)ncl-1.);
		averwcl2=(averwcl2-(wcl[clid->whichLBF]*wcl[clid->whichLBF]))/((double)ncl-1.);
		clwstats[0]=averwcl;
		clwstats[1]=averwcl2-averwcl*averwcl;
		if(clwstats[1]<=EPS*EPS){
			clwstats[1]=0.;
		}else{
			clwstats[1]=sqrt(clwstats[1]);
		}
	}else{
		clwstats[0]=wcl[0];
		clwstats[1]=0.;
	}

	free(ntotcl);
	return averw;
}
/********************************************************
*  Find the number of clusters with similar investment	*
*********************************************************/
int findSimilarInvClusters(int numcl,double targetw, double eps, double *wcl){
	int i,nscl;

	nscl=0;
	for(i=0; i<numcl; ++i){
		if((wcl[i]>=targetw-eps)&&(wcl[i]<=targetw+eps)){
			++nscl;
		}
	}

	return nscl;
}
/********************************************************
*  Calculate the difference in microbial               	*
*  composition between host @idp and its children @idk.	*
*  Store result in a vector where each element         	*
*  correspond to a different event                     	*
*********************************************************/
void storeBacDiffComp(int idp,int idk){
	int j,new;
        double mdiff;

        mdiff=0.;
        new=0;
        for(j=0; j<TYPES; ++j){
                if(bac[idp*TYPES+j]!=0.){
			mdiff+=spar->inv[j]*(bac[idp*TYPES+j]/spar->micr[idp]-bac[idk*TYPES+j]/spar->micr[idk]);
                        ++new;
                }
        }
        mdiff/=(double)new;
        offcomp->vecf[offcomp->usizef]=mdiff;
        ++offcomp->usizef;

	return;
}
/////////////////////////////////////////////////////////////////////////////////////////
/****************************************************************************************
* 				Routines for storing data on files			*
*****************************************************************************************/
/********************************************************
*   stores the average investment in the system         *
*	->Def: avInv=sum_ij(bac_ij*inv_j)/sum_ij(bacij)  *
*	where i is the host index, j is the type of     *
*   	bacteria index, bac_ij is the abundance of      *
*   	type j in host i and inv_j is the investment    *
*   	of type j.                                      *
*   	->When TYPES=2:                                  *
*   		->inv is either 0 (for neutrals)         *
*   		or 1 (for helpers)                      *
*   		->avInv is the mean frequency of         *
*   		helpers in the system                   *
*********************************************************/
void averInvestmentXt(void){
	int nh;
	double averw;
	
	nh=listh->usize;
	
	#if (NETWORK==1)//square lattice
	int nclusters;
	double clsizestats[2],clwstats[2];
	int *labels=NULL;
	int *clsize=NULL;
	double *wcl=NULL;
	ClusterMinimumID maxclsize; 
	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:t 2:<w> 3:wb 4:clsizeb 5:<wcl⁻> 6:stdwcl⁻ 7:maxclsize 8:wl⁻ 9:<clsize⁻> 10:stdclsize⁻ 11:nh/N 12:cost 13:mb 14:mh 15:rmig 16:rh\n");
		printf("#1:t 2:<w> 3:wb 4:clsizeb 5:<wcl⁻> 6:stdwcl⁻ 7:maxclsize 8:wl⁻ 9:<clsize⁻> 10:stdclsize⁻ 11:nh/N 12:cost 13:mb 14:mh 15:rmig 16:rh\n");
	}
	maxclw=malloc(sizeof(ClusterFullID));
	if (!maxclw){ perror("malloc"); exit(1);}

	labels=(int *)calloc(nh,sizeof(int));

	/*set different labels on different clusters with Hoshen-Kopelman-like algorithm using list of hosts instead of the full lattice*/
	setClusterLabelsWithList(listh,inverselisth,neighbor,LEFT,UP,labels);	

	/*fix cluster labels so they include all the integers between 0 and @numclusters-1*/
	nclusters=fixClusterLbOrder(nh,labels);
	clsize=(int *)calloc(nclusters,sizeof(int));
	/*calculating the size of each clusters and also finding the size and label of the largest cluster*/
	(void)calcClusterSizeStats(nh,nclusters,labels,clsize,clsizestats,&maxclsize);

	wcl=(double *)calloc(nclusters,sizeof(double));
	/*calculating average investment in each cluster and also finding the investment and label of the best cluster*/
	averw=findHighestInvCluster(nclusters,labels,wcl,clwstats,maxclw);//returns system average investment
	maxclw->sizeCLF=clsize[maxclw->whichLBF];//size of the best cluster
	/*storing data*/
        fprintf(gfile->file,"%f %f %f %d %f %f %d %f %f %f %f %f %f %f %d %f\n",stime->Tnow,averw,wcl[maxclw->whichLBF],clsize[maxclw->whichLBF],clwstats[0],clwstats[1],clsize[maxclsize.whichLB],wcl[maxclsize.whichLB],clsizestats[0],clsizestats[1],(double)nh/N,spar->cost,spar->mig,spar->mh,spar->rmigh,(double)spar->kh/N);
        printf("t=%f <w>=%f wb=%f clsizeb=%d <wcl⁻>=%f stdwcl⁻=%f maxclsize=%d wl=%f <clsize⁻>=%f stdclsize⁻=%f nh/N=%f cost=%f mb=%f mh=%f rmig=%d rh=%f\n",stime->Tnow,averw,wcl[maxclw->whichLBF],clsize[maxclw->whichLBF],clwstats[0],clwstats[1],clsize[maxclsize.whichLB],wcl[maxclsize.whichLB],clsizestats[0],clsizestats[1],(double)nh/N,spar->cost,spar->mig,spar->mh,spar->rmigh,(double)spar->kh/N);

	if(labels){
		free(labels);
		labels=NULL;
	}
	if(clsize){
		free(clsize);
		clsize=NULL;
	}
	if(wcl){
		free(wcl);
		wcl=NULL;
	}
	if(maxclw){
		free(maxclw);
		maxclw=NULL;
	}
	#else
	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:t 2:<w> 3:nh/N 4:cost 5:mb 6:rh\n");
		printf("#1:t 2:<w> 3:nh/N 4:cost 5:mb 6:rh\n");
	}

	averw=calcAverInv();
	fprintf(gfile->file,"%f %f %f %f %f %f\n",stime->Tnow,averw,(double)nh/N,spar->cost,spar->mig,(double)spar->kh/N);
	printf("t=%f <w>=%f nh/N=%f cost=%f mv=%f rh=%f\n",stime->Tnow,averw,(double)nh/N,spar->cost,spar->mig,(double)spar->kh/N);

	#endif
        return;
}
/**************************************
*  snapshot of host network           *
*  The colors indicate the acumulated *
*  investment of each host            *
***************************************/
void save_config(void){
        int i,j,id,idlisth,nh,cnamelen,snamelen,dnl;
        double pointsize;
	double *freqInvH=NULL;
        char *name=NULL;
	char *namedat=NULL;
	char *name_gp=NULL;
	char *nt_format=NULL;
        FILE *fconfig = NULL;
	FILE *fgp = NULL;
	
	/****creating files (if there are file with the same names, there are subscribed)********/
	
	//file name sizes
	dnl=20;
        cnamelen=gfile->fnsize+dnl;//core name length
        
	//core name structure
        name=(char *)calloc(cnamelen,sizeof(char));//base name
        sprintf(name,"snapshot_%s",gfile->fname);
	
	nt_format=(char *)calloc(dnl,sizeof(char));//format of the time term, on the name (files for animations have to be ordered: better to use a 0<t<1 format)
#if (FIG_EXT==0)
        sprintf(nt_format,"idt%f",stime->Tnow/stime->Tf);
#else
        sprintf(nt_format,"Tf%f",stime->Tnow);
#endif
	
	//file names for data and for gnuplot script
	snamelen=cnamelen+strlen(gfile->fdatapath)+strlen(nt_format)+dnl;//length of the specific file names
	namedat=(char *)calloc(snamelen,sizeof(char));//data file name
        name_gp=(char *)calloc(snamelen,sizeof(char));//name for the gnuplot data file


        sprintf(namedat,"%s%s_%s.dat",gfile->fdatapath,name,nt_format);
        sprintf(name_gp,"%s%s_%s.gp",gfile->fdatapath,name,nt_format);

        fconfig = fopen(namedat,"w");
        fgp = fopen(name_gp,"w");


	/*********Defining point size for NIntervthe gnuplot script****************************************/
        if(L>=100){
                pointsize=0.85;
        }else if(L>=50){
                pointsize=1.5;
        }else{
                pointsize=2;
        }

	/*******setting average investment per host vector and filling data file***********************/

	nh=listh->usize;
	freqInvH=(double *)calloc(nh,sizeof(double));/*indexes are the same as in @listh->vec 
						      * (real host indexes are idh=listh->vec[i] (for 0<=i<nh)) */
	calcInvFreq(freqInvH);
	
	
        for(i = 0; i < L; ++i){//linha
                for(j = 0; j < L; ++j){//coluna
                        id=j+i*L;
			idlisth=inverselisth[id];
			
			if(idlisth>=listh->usize){//@id not on the list of alive hosts
				fprintf(fconfig,"%d %d %d\n",i,j,-1);
				
			}else{//id==listh->vec[idlisth], @id is alive
				fprintf(fconfig,"%d %d %f\n",i,j,freqInvH[idlisth]);
				
			}
			
		}
                fprintf(fconfig,"\n");
	}

	/*********filling gnuplot file**********************************************/
        if(L>=100){
                pointsize=0.5;
        }else if(L>=50){
                pointsize=1.5;
        }else{
                pointsize=2;
        }
#if (FIG_EXT==0)
        fprintf(fgp,"set term png size 720,540\n");
	fprintf(fgp,"set output'%s_idt%f.png'\n",name,stime->Tnow/stime->Tf);
#else
        fprintf(fgp,"set term post eps enha color 20\n");
        fprintf(fgp,"set output'%s_Tf%f_fm.eps'\n",name,stime->Tnow);
#endif

        fprintf(fgp,"unset key\n");
        fprintf(fgp,"unset xtics\n");
        fprintf(fgp,"unset ytics\n");
        fprintf(fgp,"L=%d\n",L);
        fprintf(fgp,"set yr[0:L-1]\n");
        fprintf(fgp,"set xr[0:L-1]\n");
        fprintf(fgp,"set size square\n");
        fprintf(fgp,"set pointsize %0.2f\n",pointsize);
        fprintf(fgp,"rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b)\n");
        fprintf(fgp,"plot '%s_%s.dat' u 1:(($3>=0)?$2:1/0):(rgb(255*(1-$3),0,255*$3)) w p pt 5 lc rgb variable,\\\n",name,nt_format);
        fprintf(fgp,"'%s_%s.dat' u 1:(($3<0)?$2:1/0):(rgb(0,0,0)) w p pt 5 lc rgb variable\n",name,nt_format);


	/***freeing allocated memory and closing files******/
        free(freqInvH);
	
	free(nt_format);
        free(name);
        free(name_gp);
        free(namedat);
	
	fclose(fconfig);
        fclose(fgp);

        return;
}
/*****************************************************
 * store EvMeasures *meas;investment distribution among hosts and     *
 * create a gnuplot script to create graphics        *
 *****************************************************/
void invDistXt(void){
	int i,nh,nbins=100,namelen,dnl;
	double binsize,*hist_inv,*freqInvH;
        char *name,*namedat,*name_gp,*nt_format;
        FILE *fhist,*fgp;

	/****creating files (if there are file with the same names, there are subscribed)********/
	dnl=100;
	namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;
        nt_format=(char *)calloc(dnl,sizeof(char));//base name
        name=(char *)calloc(namelen,sizeof(char));//format of the time term, on the name (files for animations have to be ordered: better to use a 0<t<1 format)
        namedat=(char *)calloc((namelen+dnl),sizeof(char));//data file name
        name_gp=(char *)calloc((namelen+dnl),sizeof(char));//name for the gnuplot data file
        
	sprintf(name,"invDistXt_%s",gfile->fname);

#if (FIG_EXT==0)
        sprintf(nt_format,"idt%f",(double)stime->Tnow/stime->Tf);
#else
        sprintf(nt_format,"Tf%f",stime->Tnow);
#endif
        sprintf(namedat,"%s%s_%s.dat",gfile->fdatapath,name,nt_format);
        sprintf(name_gp,"%s%s_%s.gp",gfile->fdatapath,name,nt_format);
        
	fhist = fopen(namedat,"w");
        fgp = fopen(name_gp,"w");
	
	/*******setting histogram and investment frequency per host vectors***********************/
	
	hist_inv=(double *)calloc(nbins,sizeof(double));
	memset(hist_inv,0.,sizeof(double)*nbins);
	binsize=1./(double)nbins;
	
	nh=listh->usize;
	freqInvH=(double *)calloc(nh,sizeof(double));/*indexes are the same as in @listh->vec 
						      * (real host indexes are idh=listh->vec[i] (for 0<=i<nh)) */
	calcInvDist(binsize,hist_inv,freqInvH);

	/*******************filling data file*****************************/
	
	if(stime->Tnow==0.)fprintf(fhist,"#1:average inv. in hosts 2:frac. of hosts 3:numsteps\n");
	for(i=0; i<nbins; ++i){
		fprintf(fhist,"%f %f %f %f\n",(double)i*binsize,(double)hist_inv[i]/nh,freqInvH[0],stime->Tnow);
	}
	
	/*******************filling gnuplot file*****************************/

#if (FIG_EXT==0)
       fprintf(fgp,"set term png size 720,540\n");
       fprintf(fgp,"set output'%s_idt%f.png'\n",name,(double)stime->Tnow/stime->Tf);
#else
        fprintf(fgp,"set term post eps enha color 20\n");
        fprintf(fgp,"set output'%s_Tf%f_fm.eps'\n",name,stime->Tnow);
#endif

	fprintf(fgp,"set title'{/=15 time steps=%f}'\n",stime->Tnow);
	fprintf(fgp,"set xr[0:1]\n");
	fprintf(fgp,"set yr[0:1]\n");
	fprintf(fgp,"set xlabel'{/=25 av. inv. in host}'\n");
	fprintf(fgp,"set ylabel'{/=25 frac. of hosts}'\n");
	fprintf(fgp,"plot '%s_%s.dat' u 1:2 w l lw 3 lc 1\n",name,nt_format);

	/**freeing allocated memory and closing files**/

	free(freqInvH);
	free(hist_inv);
	free(name);
	free(namedat);
	free(name_gp);
	free(nt_format);
	fclose(fhist);
	fclose(fgp);
	return;
}
/************************************************************************
* 	Stores the following cluster related measures: 			*
* 	->Highest cluster average investment (=best clluster)		*
* 	->Cluster size of the best cluster				*
* 	->Investment of the largest cluster				*
* 	->Size of the largest cluster					*
* 	->number of clusters						*
* 	->average cluster size and related standart deviation		*
* 	->System average investment					*
* 	->Number of hosts						*		
*************************************************************************/
void clustersXt(void){
	int i,nh,nclusters;
	double averw,clstats[2],clwstats[2],meanwcl,meanwcl2,stdwcl;
	int *labels=NULL;
	int *clsize=NULL;
	double *wcl=NULL;
	ClusterMinimumID maxclsize; 
	
	nh=listh->usize;

	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:time 2:<w> of the best cluster 3:cluster size of the best cluster 4:<w> of the largest cluster 5:size of the largest cluster 6:#ofclusters 7:mean cluster size 8:std(clsize) 9:<w> 10:<wcl> 11:std(wcl) 12:nh\n");
	}
	maxclw=malloc(sizeof(ClusterFullID));
	if (!maxclw){ perror("malloc"); exit(1);}

	labels=(int *)calloc(nh,sizeof(int));

	/*set different labels on different clusters with Hoshen-Kopelman-like algorithm using list of hosts instead of the full lattice*/
	setClusterLabelsWithList(listh,inverselisth,neighbor,LEFT,UP,labels);	

	/*fix cluster labels so they include all the integers between 0 and @numclusters-1*/
	nclusters=fixClusterLbOrder(nh,labels);
	clsize=(int *)calloc(nclusters,sizeof(int));
	/*calculating the size of each clusters and also finding the size and label of the largest cluster*/
	calcClusterSizeStats(nh,nclusters,labels,clsize,clstats,&maxclsize);

	wcl=(double *)calloc(nclusters,sizeof(double));
	/*calculating average investment in each cluster and also finding the investment and label of the best cluster*/
	averw=findHighestInvCluster(nclusters,labels,wcl,clwstats,maxclw);//returns system average investment
	maxclw->sizeCLF=clsize[maxclw->whichLBF];//size of the best cluster

	/*calculating mean cluster average investment and the related standart deviation*/	
	meanwcl=0.;
	meanwcl2=0.;
	for(i=0; i<nclusters; ++i){
		meanwcl+=wcl[i];
		meanwcl2+=wcl[i]*wcl[i];
	}
	meanwcl/=(double)nclusters;
	meanwcl2/=(double)nclusters;
	if(meanwcl2>0.){
		stdwcl=meanwcl2-meanwcl*meanwcl;
		stdwcl=sqrt(stdwcl);
	}else{
		stdwcl=0.;
	}
	/*storing data on file*/
        fprintf(gfile->file,"%f %f %d %f %d %d %f %f %f %f %f %d\n",stime->Tnow,maxclw->stateCL,maxclw->sizeCLF,wcl[maxclsize.whichLB],maxclsize.sizeCL,nclusters,clstats[0],clstats[1],averw,meanwcl,stdwcl,nh);
        printf("t=%f bestClInv=%f BestClsize=%d LargestClInv=%f LargestClsize=%d ncl=%d <clsize>=%f std(clsize)=%f <w>=%f <wcl>=%f std(wcl)=%f nh=%d\n",stime->Tnow,maxclw->stateCL,maxclw->sizeCLF,wcl[maxclsize.whichLB],maxclsize.sizeCL,nclusters,clstats[0],clstats[1],averw,meanwcl,stdwcl,nh);

	/****freeing allocated memory*****/
	free(labels);
	free(wcl);
	free(clsize);
	if(maxclw){
		free(maxclw);
		maxclw=NULL;
	}
	return;
}
/********************************************************
*       Stores transient investment distribution among 	*
*       clusters: 					*
*       -> time 					*
*       -> all possible average investments, <w>	*
*       ->fraction of clusters in each 			*
*         bin (int)(<w>/binsize)			*
*********************************************************/
void clustersInvDistXt(void){
	int nh,nclusters;
	int i,id,idmax,nbins=100;
	double averw,clwstats[2];
	double binsize;
	int *hist=NULL;
	int *labels=NULL;
	double *wcl=NULL;
	
	binsize=1./(double)nbins;//investment bin size

	nh=listh->usize;

	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:time 2:<w> 3:frac. of clusters with average investment=<w> |Binsize=%f\n",binsize);
	}

	labels=(int *)calloc(nh,sizeof(int));
	maxclw=malloc(sizeof(ClusterFullID));
	if (!maxclw){ perror("malloc"); exit(1);}

	/*set different labels on different clusters with Hoshen-Kopelman-like algorithm using list of hosts instead of the full lattice*/
	setClusterLabelsWithList(listh,inverselisth,neighbor,LEFT,UP,labels);	

	/*fix cluster labels so they include all the integers between 0 and @numclusters-1*/
	nclusters=fixClusterLbOrder(nh,labels);

	wcl=(double *)calloc(nclusters,sizeof(double));
	/*calculating average investment in each cluster and also finding the investment and label of the best cluster*/
	averw=findHighestInvCluster(nclusters,labels,wcl,clwstats,maxclw);
	/*Histogram: calculating the fraction of clusters for each investment bin*/
	idmax=nbins-1;
	hist=(int *)calloc(nbins,sizeof(int));
	memset(hist,0,sizeof(int)*nbins);
	for(i=0; i<nclusters;++i){
		id=(int)(wcl[i]/binsize);
		if(id>idmax)id=idmax;
		++hist[id];
	}
	/*storing data*/
	for(i=0; i<nbins; ++i){
		fprintf(gfile->file,"%f %f %f %f %d\n",stime->Tnow,(double)i/nbins,(double)hist[i]/nclusters,averw,nh);
		printf("t=%f <w>=%f cl_frac=%f <w>sys=%f nh=%d\n",stime->Tnow,(double)i/nbins,(double)hist[i]/nclusters,averw,nh);
	}
	fprintf(gfile->file,"%f %f %f %f %d\n\n",stime->Tnow,1.,(double)hist[nbins-1]/nclusters,averw,nh);

	/****freeing allocated memory*****/
	free(hist);
	free(labels);
	free(wcl);
	if(maxclw){
		free(maxclw);
		maxclw=NULL;
	}
	return;
}
/********************************************************
 * Measures related to the stability of the cluster	*
 * with the highest average investment.			*
 * 							*
 * Stores: time of emergence of a cluster with high	*
 * investment (>=@BESTWtr) and that is also the 	*
 * cluster with highest average investment; time 	*
 * of where this it is no longer the best cluster 	*
*********************************************************/
void bestClusterStabilityMeasures(void){
	int i,nh,nclusters,nblist,numsimilarw;
	double averw, clstats[2],clwstats[2],eps;
	int *clsize=NULL;
	double *wcl=NULL;
	ClusterMinimumID maxclsize; 
	
	nh=listh->usize;

	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:t 2:wcl_best 3:clsize 4:<w> 5:#of clusters with wcl~wcl_best 6:#of clusters  participants: \n");
	}

	/*set different labels on different clusters with Hoshen-Kopelman-like algorithm using list of hosts instead of the full lattice*/
	setClusterLabelsWithList(listh,inverselisth,neighbor,LEFT,UP,lb);	

	/*fix cluster labels so they include all the integers between 0 and @numclusters-1*/
	nclusters=fixClusterLbOrder(nh,lb);
	clsize=(int *)calloc(nclusters,sizeof(int));
	/*calculating the size of each clusters and also finding the size and label of the largest cluster*/
	calcClusterSizeStats(nh,nclusters,lb,clsize,clstats,&maxclsize);

	wcl=(double *)calloc(nclusters,sizeof(double));
	/*calculating average investment in each cluster and also finding the investment and label of the best cluster*/
	averw=findHighestInvCluster(nclusters,lb,wcl,clwstats,maxclw);//returns system average investment
	maxclw->sizeCLF=clsize[maxclw->whichLBF];//size of the best cluster

	eps=spar->inv[1];//=dw between the investments of a type of bacteria j and another of type j+1 or j-1: dw=0.01 for T=101
	numsimilarw=findSimilarInvClusters(nclusters,maxclw->stateCL,eps,wcl);
	/**************/
	nblist=0;

	for(i=0; i<nh; ++i){
		if(lb[i]==maxclw->whichLBF){
			bestwlisth->vec[nblist]=i;
			++nblist;
		}
		bestwlisth->usize=nblist;
	}
	fprintf(gfile->file,"%f %f %d %f %d %d ",stime->Tnow,wcl[maxclw->whichLBF],clsize[maxclw->whichLBF],averw,numsimilarw,nclusters);
	printf("t=%f wclb=%f clsizeb=%d <w>=%f nswcl=%d numcl=%d partic.: ",stime->Tnow,wcl[maxclw->whichLBF],clsize[maxclw->whichLBF],averw,numsimilarw,nclusters);
	for(i=0; i<nblist; ++i){
		fprintf(gfile->file,"%d ",bestwlisth->vec[i]);
		printf("%d ",bestwlisth->vec[i]);
	}
	fprintf(gfile->file,"\n");
	printf("\n");

	/****freeing allocated memory*****/
	free(wcl);
	free(clsize);
	return;
}
/***************************************************
*  stores the number of host events per            * 
*  microbial time steps                            *
****************************************************/
void numHostEventsPerDtXt(void){

	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:time 2:#of host events per Dt_ref 3:#of host births per Dt_ref 4:#of host deaths per Dt_ref 5:#of hosts\n");
	}
	fprintf(gfile->file,"%f %d %d %d %d\n",stime->Tnow,meas->numb+meas->numd,meas->numb,meas->numd,listh->usize);
	printf("%f %d %d %d %d\n",stime->Tnow,meas->numb+meas->numd,meas->numb,meas->numd,listh->usize);

	return;
}
/****************************************************************
* 	store spatial correlation for lattice ocupation 	*
* 	as a function of time and distance 			*
*****************************************************************/
void spatialCorrXt(void){
	int i,j,idh,d,dmax,idx,idy,nh;
	double corrh,averinv,tot_micr;

	if(stime->Tnow==0.){
		if(gfile->file==NULL){
			printf("You are trying to write in a file that doesn't exist.\n");
			exit(1);
		}else{
			fprintf(gfile->file,"#1:time 2:Distance 3:corrh 4:<w> 5:ntot 6:nh/N 7:nh\n");
			printf("#1:time 2:Distance 3:corrh 4:<w> 5:ntot 6:nh/N 7:nh\n");
		}
	}

	nh=listh->usize;

	averinv=calcAverInv();
	tot_micr=0.;
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		tot_micr+=spar->micr[idh];
	}

	dmax=L/2;
	for(d=1; d<=dmax; ++d){
		corrh=0.;
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			idx=idh;
			idy=idh;
			for(j=0; j<d; ++j){
				idx=neighbor[idx][RIGHT];
				idy=neighbor[idy][DOWN];
			}
			//host[i]*host[j]=1, if host[i]=host[j]=1, or 0 otherwise
			corrh+=(double)(host[idh]*host[idx]+host[idh]*host[idy]);

		}
		corrh/=(2.*nh);

		fprintf(gfile->file,"%f %d %f %f %f %f %d\n",stime->Tnow,d,corrh,averinv,tot_micr,(double)listh->usize/N,listh->usize);
		printf("t=%f d=%d corrh=%f <w>=%f ntot=%f rho_o=%f nh=%d\n",stime->Tnow,d,corrh,averinv,tot_micr,(double)listh->usize/N,listh->usize);
	}
	fprintf(gfile->file,"\n");
	printf("\n");

	return;
}
/********************************************************
*  	Measure average difference of hosts parents	* 
*  	and their childrens microbial composition 	*
*  	over time (to compare different vertical 	*
*  	transmission approachs)      			*
*********************************************************/
void difMicrCompXt(void){
	int numb;
	int nh=listh->usize;
	double averinv;
	double stats[2];


	if(stime->Tnow==0.){
		#if (OFFCOMP==0)
		fprintf(gfile->file,"#1:time 2:average dif. between parent-offspring freqInv 3:stardart deviation 4:average invest. 5:sample (#of repr.) 6:nh 7:nh/N\n");
		#else
		//<w[offspring]/micr[offcalcAverInvspring]> is the mean of the acumulated investment divided by the amount of microbes received by newborns (average of last @SAMPLE reproductions)
		fprintf(gfile->file,"#1:time 2:<w[offpring]/micr[offspring]>  3:stardart deviation 4:average invest. 5:sample (#of repr.) 6:nh 7:nh/N\n");
		#endif
	}

	numb =offcomp->usizef;
	if(numb>=SAMPLE){
		calcRMSError(offcomp->vecf,0,offcomp->usizef-1,stats);
	
		averinv=calcAverInv();
		
		fprintf(gfile->file,"%f %f %f %f %d %d %f\n",stime->Tnow,stats[0],stats[1],averinv,numb,nh,(double)nh/N);
		printf("t=%f mean=%f std=%f avinv=%f nb=%d nh=%d nh/N=%f\n",stime->Tnow,stats[0],stats[1],averinv,numb,nh,(double)nh/N);

		//reset	
		memset(offcomp->vecf,0.,sizeof(double)*offcomp->sizef);
		offcomp->usizef=0;
	}
	
	return;
}
/************************************************************************
*  Store in @SAMPLE files the average investment over time		*
*  for multiple costs (@sample files for each cost).			*
*  The routine changes the values of the cost, but the measuring and 	*
*  storing is made by @averInvestmentXt(), which is being called in the 	*
*  the evolution routine @callSysDynamics().				*	
*************************************************************************/
void averInvXtMultipleCosts(Event *event,Event *mevent){
        int i;
	double dc, costmax;

        stime->saveT=0.;
        stime->Tf=TF;
        
	spar->cost=0.01;
	dc=0.01;
	costmax=0.2;

	while(spar->cost<=costmax){
		allocateMemTM();
        #if (NETWORK==0)
		if(spar->cost<=0.1){
			stime->Tf=60000;
		}else{
			stime->Tf=120000;
		}
	#else
		if(spar->cost<=0.05){
			stime->Tf=20000;
		}else if(spar->cost<0.1){
			stime->Tf=80000;
		}else if(spar->cost<=0.15){
			stime->Tf=100000;
		}else{
			stime->Tf=150000;
		}
	#endif
		for(i=0; i<SAMPLE; ++i){
			openFiles();
			callSysDynamics(event,mevent);
			closeFiles();
			stime->Tnow=0.;
			stime->saveT=0.;
			setCI();
		}
		spar->cost+=dc;
		freeMemTM();
	}


	return;
}
/*******steady state measures***********************************************************/

/****************************************************************
*  Store the steady state average investment for different  	*
*  costs and probability values of long range migration, plr	*	
*****************************************************************/
void costXplrXw(Event *event,Event *mevent){
	int i,nh;
        int ok=0,namelen,dnl;
	double dplr,dc,cmax,cmin,plrmax,stats[2];
	int npar=2,mul,nf;
	double expo,param[npar];
	unsigned long id;
        char nparam[npar][10];
	
	/******seting file*************************************/ 
        //generic file struct
	gfile=malloc(sizeof(GenFile));
	if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=300;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");
	
	//file name components
        param[0]=spar->mu;
        param[1]=spar->mig;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
			expo=floor(log10(param[i]));
        		mul=ceil(param[i]/pow(10.,expo));
			nf=0;
			while(ceil(param[i]/pow(10.,expo))!=(int)(param[i]/pow(10.,expo))){
				--expo;
				mul=ceil(param[i]/pow(10.,expo));
				++nf;
			}
                        sprintf(nparam[i],"%de%d",mul,(int)expo);
                }
        }

        sprintf(gfile->fname,"SL_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_mh%0.1f_rmh%d_plr%0.2f_MT%d",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->mh,spar->rmigh,spar->plr,MIGRATION_TYPE);
        dnl=200;
	namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;
	char *name=(char *)calloc(namelen,sizeof(char));

        id = (unsigned long)time(NULL);
        
        while(ok==0){
		sprintf(name,"%scostXplrXw_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
               
	sprintf(name,"%scostXplrXw_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
	
	/****Allocate memory**********************************************/

	avinv=malloc(sizeof(DynVec));
	if (!avinv) { perror("malloc"); exit(1);}
	avinv->sizef=40000;
	avinv->usizef=0;
	avinv->vecf=(double *)calloc(avinv->sizef,sizeof(double));

	cmax=0.2;
	cmin=1e-2;
	spar->cost=cmin;
	dc=cmin;
	dplr=0.1;
	spar->plr=0.;
	plrmax=1.;

	while(spar->cost<=cmax){
		if(spar->cost<=0.05){
			stime->Tf=20000;
		}else if(spar->cost<=0.07){
			stime->Tf=30000;
		}else if(spar->cost<=0.09){
			stime->Tf=40000;
		}else if(spar->cost<=0.11){
			stime->Tf=80000;
		}else if(spar->cost<=0.15){
			stime->Tf=100000;
		}else{
			stime->Tf=150000;
		}
		
		stime->timewindow=ceil(0.25*stime->Tf);
		stime->transtime=stime->Tf-stime->timewindow;
	
		spar->plr=0.;
		while(spar->plr<=plrmax){
			setCI();
			stime->Tnow=0.;
			stime->saveT=stime->transtime;
			callSysDynamics(event,mevent);
			
			nh=listh->usize;

			calcRMSError(avinv->vecf,0,avinv->usizef,stats);
		
			fprintf(gfile->file,"%f %f %f %f %d %f\n",spar->cost,spar->plr,stats[0],stats[1],nh,stime->Tf);
			fflush(gfile->file);
			printf("cost=%f plr=%f <w>=%f sdtw=%f nh=%d Tf=%f\n",spar->cost,spar->plr,stats[0],stats[1],nh,stime->Tf);
		
			avinv->usizef=0;

			spar->plr+=dplr;
		}
		fprintf(gfile->file,"\n");
		printf("\n");

		spar->cost+=dc;
	}
	
	if(gfile->file){
		fclose(gfile->file);
		gfile->file=NULL;
	}
	if(gfile->fname){
		free(gfile->fname);
		gfile->fname=NULL;
	}
	if(gfile->fdatapath){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile){
		free(gfile);
		gfile=NULL;
	}

	free(name);
	return;
}
/****************************************************************
*  Store in @SAMPLE files the average investment		*
*  in the system as a function of the system carrying dilution	*	
*****************************************************************/
void averInvXrh(Event *event,Event *mevent){
	int i,idh,nh,steady;
        int ok=0,namelen,dnl;
	double rh,drh=0.05,eps=0.05,tot_micr,stats[2];
	double twind=stime->timewindow;	
	int npar=2,mul,nf;
	double expo,param[npar];
	unsigned long id;
        char nparam[npar][10];
	
       
	/******seting file*************************************/ 
        //generic file struct
	gfile=malloc(sizeof(GenFile));
	if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=400;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");
	
	//file name components
        param[0]=spar->mu;
        param[1]=spar->mig;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
			expo=floor(log10(param[i]));
        		mul=ceil(param[i]/pow(10.,expo));
			nf=0;
			while(ceil(param[i]/pow(10.,expo))!=(int)(param[i]/pow(10.,expo))){
				--expo;
				mul=ceil(param[i]/pow(10.,expo));
				++nf;
			}
                        sprintf(nparam[i],"%de%d",mul,(int)expo);
                }
        }

	#if(NETWORK==0)
        sprintf(gfile->fname,"CG_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1]);
	#else
        sprintf(gfile->fname,"SL_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_mh%0.1f_rmh%d_plr%0.2f_MT%d",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->mh,spar->rmigh,spar->plr,MIGRATION_TYPE);
	#endif
        dnl=200;
	namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;
	char *name=(char *)calloc(namelen,sizeof(char));

        id = (unsigned long)time(NULL);
        
        while(ok==0){
                sprintf(name,"%savInvXrh_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
               
	sprintf(name,"%savInvXrh_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
	
	/****Allocate memory**********************************************/

	avinv=malloc(sizeof(DynVec));
	if (!avinv) { perror("malloc"); exit(1);}
	avinv->sizef=stime->timewindow/stime->tinterval+1;
	avinv->usizef=0;
	avinv->vecf=(double *)calloc(avinv->sizef,sizeof(double));

	/****Dynamics******************/

	rh=drh;
	spar->kh=N*rh;

	while(rh<=0.75){
		setCI();

		stime->saveT=stime->transtime;
		stime->Tf=stime->saveT+twind;
		steady=0;
		do{
			callSysDynamics(event,mevent);

			calcRMSError(avinv->vecf,0,avinv->usizef,stats);//test if system reached equilibrium
			if(stats[1]<eps){
				steady=1;
			}else{
				stime->Tf+=twind;
				stime->saveT=stime->Tnow;
				avinv->usizef=0.;
			}
			printf("steady=%d time=%f <avinv>=%f std=%f\n",steady,stime->Tnow,stats[0],stats[1]);
		}while(steady==0);

		nh=listh->usize;
		tot_micr=0.;
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			tot_micr+=spar->micr[idh];
		}
		
		fprintf(gfile->file,"%f %f %f %d %f %f\n",rh,stats[0],stats[1],nh,stime->Tf,tot_micr);
		fflush(gfile->file);
		printf("%f %f %f %d %f %f\n",rh,stats[0],stats[1],nh,stime->Tf,tot_micr);
		
		avinv->usizef=0;

		rh+=drh;
		spar->kh=N*rh;
	}
		
	/***freeing memory*****/
	free(avinv->vecf);
	free(avinv);

	if(gfile->file){
		fclose(gfile->file);
		gfile->file=NULL;
	}
	if(gfile->fname){
		free(gfile->fname);
		gfile->fname=NULL;
	}
	if(gfile->fdatapath){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile){
		free(gfile);
		gfile=NULL;
	}

	free(name);
	return;
}
/************************************************************************
*  Store in @SAMPLE files the average investment			*
*  in the system as a function of the system number of microbial	*
*  generations per host generations, Gh					*	
*************************************************************************/
void averInvXgh(Event *event,Event *mevent){
	int i,idh,nh,steady,gh_max,dgh=5;
        int ok=0,namelen,dnl;
	double eps=0.05,tot_micr,stats[2];
	double twind=stime->timewindow;	
	int npar=2,mul,nf;
	double expo,param[npar];
	unsigned long id;
        char nparam[npar][10];
	
       
	/******seting file*************************************/ 
        //generic file struct
	gfile=malloc(sizeof(GenFile));
	if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=400;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");
	
	//file name components
        param[0]=spar->mu;
        param[1]=spar->mig;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
			expo=floor(log10(param[i]));
        		mul=ceil(param[i]/pow(10.,expo));
			nf=0;
			while(ceil(param[i]/pow(10.,expo))!=(int)(param[i]/pow(10.,expo))){
				--expo;
				mul=ceil(param[i]/pow(10.,expo));
				++nf;
			}
                        sprintf(nparam[i],"%de%d",mul,(int)expo);
                }
        }

	#if(NETWORK==0)
        sprintf(gfile->fname,"CG_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1]);
	#else
        sprintf(gfile->fname,"SL_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_mh%0.1f_rmh%d_plr%0.2f_MT%d",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->mh,spar->rmigh,spar->plr,MIGRATION_TYPE);
	#endif
        dnl=200;
	namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;
	char *name=(char *)calloc(namelen,sizeof(char));

        id = (unsigned long)time(NULL);
        
        while(ok==0){
                sprintf(name,"%savInvXgh_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
               
	sprintf(name,"%savInvXgh_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
	
	/****Allocate memory**********************************************/

	avinv=malloc(sizeof(DynVec));
	if (!avinv) { perror("malloc"); exit(1);}
	avinv->sizef=stime->timewindow/stime->tinterval+1;
	avinv->usizef=0;
	avinv->vecf=(double *)calloc(avinv->sizef,sizeof(double));

	/****Dynamics******************/

	spar->gh=Gh;
	gh_max=Gh+200;

	while(spar->gh<=gh_max){
		setCI();

		stime->saveT=stime->transtime;
		stime->Tf=stime->saveT+twind;
		steady=0;
		do{
			callSysDynamics(event,mevent);

			calcRMSError(avinv->vecf,0,avinv->usizef,stats);//test if system reached equilibrium
			if(stats[1]<eps){
				steady=1;
			}else{
				stime->Tf+=twind;
				stime->saveT=stime->Tnow;
				avinv->usizef=0.;
			}
			printf("steady=%d time=%f <avinv>=%f std=%f\n",steady,stime->Tnow,stats[0],stats[1]);
		}while(steady==0);

		nh=listh->usize;
		tot_micr=0.;
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			tot_micr+=spar->micr[idh];
		}
		
		fprintf(gfile->file,"%d %f %f %d %f %f\n",spar->gh,stats[0],stats[1],nh,stime->Tf,tot_micr);
		fflush(gfile->file);
		printf("%d %f %f %d %f %f\n",spar->gh,stats[0],stats[1],nh,stime->Tf,tot_micr);
		
		avinv->usizef=0;

		spar->gh+=dgh;
	}
		
	/***freeing memory*****/
	free(avinv->vecf);
	free(avinv);

	if(gfile->file){
		fclose(gfile->file);
		gfile->file=NULL;
	}
	if(gfile->fname){
		free(gfile->fname);
		gfile->fname=NULL;
	}
	if(gfile->fdatapath){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile){
		free(gfile);
		gfile=NULL;
	}
	free(name);
	return;
}
/************************************************************************
*  Store in @SAMPLE files the average investment			*
*  in the system as a function of the microbial migration rate		*	
*************************************************************************/
void averInvXmb(Event *event,Event *mevent,double mbmin,double mbmax){
	int i,idh,nh,steady;
        int ok=0,namelen,dnl;
	double dmb,eps=0.05,tot_micr,stats[2];
	double twind=stime->timewindow;	
	int mul,nf;
	double expo,param;
	unsigned long id;
        char nparam[10];
	
       
	/******seting file*************************************/ 
        //generic file struct
	gfile=malloc(sizeof(GenFile));
	if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=300;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");
	
	//file name components
        param=spar->mu;

	if(param==0.){
		sprintf(nparam,"0");
	}else{
		expo=floor(log10(param));
		mul=ceil(param/pow(10.,expo));
		nf=0;
		while(ceil(param/pow(10.,expo))!=(int)(param/pow(10.,expo))){
			--expo;
			mul=ceil(param/pow(10.,expo));
			++nf;
		}		
		sprintf(nparam,"%de%d",mul,(int)expo);
	}

	#if(NETWORK==0)
        sprintf(gfile->fname,"CG_L%d_Ty%d_Kh%d_cost%0.2f_mu%s",L,TYPES,spar->kh,spar->cost,nparam);
	#else
        sprintf(gfile->fname,"SL_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mh%0.1f_rmh%d_plr%0.2f_MT%d",L,TYPES,spar->kh,spar->cost,nparam,spar->mh,spar->rmigh,spar->plr,MIGRATION_TYPE);
	#endif
        dnl=200;
	namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;
	char *name=(char *)calloc(namelen,sizeof(char));

        id = (unsigned long)time(NULL);
        
        while(ok==0){
                sprintf(name,"%savInvXmb_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
        
	sprintf(name,"%savInvXmb_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
	
	/****Allocate memory**********************************************/

	avinv=malloc(sizeof(DynVec));
	if (!avinv) { perror("malloc"); exit(1);}
	avinv->sizef=stime->timewindow/stime->tinterval+1;
	avinv->usizef=0;
	avinv->vecf=(double *)calloc(avinv->sizef,sizeof(double));

	/****Dynamics******************/

	spar->mig=mbmin;

	while(spar->mig<=mbmax){
		setCI();

		stime->saveT=stime->transtime;
		stime->Tf=stime->saveT+twind;
		steady=0;
		do{
			callSysDynamics(event,mevent);

			calcRMSError(avinv->vecf,0,avinv->usizef,stats);//test if system reached equilibrium
			if(stats[1]<eps){
				steady=1;
			}else{
				stime->Tf+=twind;
				stime->saveT=stime->Tnow;
				avinv->usizef=0.;
			}
			printf("steady=%d time=%f <avinv>=%f std=%f\n",steady,stime->Tnow,stats[0],stats[1]);
		}while(steady==0);

		nh=listh->usize;
		tot_micr=0.;
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			tot_micr+=spar->micr[idh];
		}
		
		fprintf(gfile->file,"%f %f %f %d %f %f\n",spar->mig,stats[0],stats[1],nh,stime->Tf,tot_micr);
		fflush(gfile->file);
		printf("mb=%f <w>=%f sdt=%f nh=%d Tf=%f micr=%f\n",spar->mig,stats[0],stats[1],nh,stime->Tf,tot_micr);
		
		avinv->usizef=0;
		
		#if (NETWORK==0)
			if(spar->mig<1e-5){
				dmb=5e-6;
			}else if(spar->mig<0.002){
				dmb=1e-4;
			}else{
				dmb=5e-3;
			}
		#else
			if(spar->mig<1e-5){
				dmb=5e-6;
			}else if(spar->mig<4e-4){
				dmb=1e-5;
			}else{
				dmb=1e-5;
			}if(spar->mig<1e-3){
				dmb=1e-4;
			}else{
				dmb=5e-3;
			}
		#endif
		spar->mig+=dmb;
	}

	/***freeing memory*****/
	free(avinv->vecf);
	free(avinv);

	if(gfile->file){
		fclose(gfile->file);
		gfile->file=NULL;
	}
	if(gfile->fname){
		free(gfile->fname);
		gfile->fname=NULL;
	}
	if(gfile->fdatapath){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile){
		free(gfile);
		gfile=NULL;
	}
	free(name);
	return;
}
/************************************************************************
*  Store in @SAMPLE files the average investment			*
*  in the system as a function of the cost factor payed by helpfull	* 
*  bacteria								*	
*************************************************************************/
void averInvXcost(Event *event,Event *mevent){
	int i,idh,nh,steady;
        int ok=0,namelen,dnl;
	double eps=0.05,tot_micr,stats[2],cost_max,dc;
	double twind=stime->timewindow;	
	int npar=2,mul,nf;
	double param[npar],expo;
	unsigned long id;
        char nparam[npar][10];
	
       
	/******seting file*************************************/ 
        //generic file struct
	gfile=malloc(sizeof(GenFile));
	if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=300;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");
	
	//file name components
        param[0]=spar->mu;
        param[1]=spar->mig;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
			expo=floor(log10(param[i]));
        		mul=ceil(param[i]/pow(10.,expo));
			nf=0;
			while(ceil(param[i]/pow(10.,expo))!=(int)(param[i]/pow(10.,expo))){
				--expo;
				mul=ceil(param[i]/pow(10.,expo));
				++nf;
			}
                        sprintf(nparam[i],"%de%d",mul,(int)expo);
                }
        }

	#if (NETWORK==0)
	sprintf(gfile->fname,"CG_L%d_Ty%d_Kh%d_mu%s_mb%s",L,TYPES,spar->kh,nparam[0],nparam[1],spar->mh,spar->rmigh);
	#else
	sprintf(gfile->fname,"SL_L%d_Ty%d_Kh%d_mu%s_mb%s_mh%0.1f_rmh%d_plr%0.2f_MT%d",L,TYPES,spar->kh,nparam[0],nparam[1],spar->mh,spar->rmigh,spar->plr,MIGRATION_TYPE);
	#endif
        dnl=200;
	namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;
	char *name=(char *)calloc(namelen,sizeof(char));

        id = (unsigned long)time(NULL);
        
        while(ok==0){
                sprintf(name,"%savInvXcost_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
               
	sprintf(name,"%savInvXcost_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
	
	/****Allocate memory**********************************************/

	avinv=malloc(sizeof(DynVec));
	if (!avinv) { perror("malloc"); exit(1);}
	avinv->sizef=stime->timewindow/stime->tinterval+1;
	avinv->usizef=0;
	avinv->vecf=(double *)calloc(avinv->sizef,sizeof(double));

	/****Dynamics******************/

	spar->cost=1e-2;
	cost_max=0.21;
	dc=0.01;


	while(spar->cost<=cost_max){
		setCI();

		stime->saveT=stime->transtime;
		stime->Tf=stime->saveT+twind;
		steady=0;
		do{
			callSysDynamics(event,mevent);

			calcRMSError(avinv->vecf,0,avinv->usizef,stats);//test if system reached equilibrium
			if(stats[1]<eps){
				steady=1;
			}else{
				stime->Tf+=twind;
				stime->saveT=stime->Tnow;
				avinv->usizef=0.;
			}
			printf("steady=%d time=%f <avinv>=%f std=%f\n",steady,stime->Tnow,stats[0],stats[1]);
		}while(steady==0);

		nh=listh->usize;
		tot_micr=0.;
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			tot_micr+=spar->micr[idh];
		}
		
		fprintf(gfile->file,"%f %f %f %d %f %f\n",spar->cost,stats[0],stats[1],nh,stime->Tf,tot_micr);
		fflush(gfile->file);
		printf("%f %f %f %d %f %f\n",spar->cost,stats[0],stats[1],nh,stime->Tf,tot_micr);
		
		avinv->usizef=0;

		spar->cost+=dc;
	}
		
	/***freeing memory*****/
	free(avinv->vecf);
	free(avinv);

	if(gfile->file){
		fclose(gfile->file);
		gfile->file=NULL;
	}
	if(gfile->fname){
		free(gfile->fname);
		gfile->fname=NULL;
	}
	if(gfile->fdatapath){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile){
		free(gfile);
		gfile=NULL;
	}
	free(name);
	return;
}
/************************************************************************
*  Store in @SAMPLE files the average investment			*
*  in the system as a function of the fraction of long-range host 	*
*  migration 								*	
*************************************************************************/
void averInvXplr(Event *event,Event *mevent){
	int i,idh,nh,steady;
        int ok=0,namelen,dnl;
	double eps=0.05,tot_micr,stats[2],plrmax,dp;
	double twind=stime->timewindow;	
	int npar=2,mul,nf;
	double param[npar],expo;
	unsigned long id;
        char nparam[npar][10];
	
       
	/******seting file*************************************/ 
        //generic file struct
	gfile=malloc(sizeof(GenFile));
	if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=300;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");
	
	//file name components
        param[0]=spar->mu;
        param[1]=spar->mig;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
			expo=floor(log10(param[i]));
        		mul=ceil(param[i]/pow(10.,expo));
			nf=0;
			while(ceil(param[i]/pow(10.,expo))!=(int)(param[i]/pow(10.,expo))){
				--expo;
				mul=ceil(param[i]/pow(10.,expo));
				++nf;
			}
                        sprintf(nparam[i],"%de%d",mul,(int)expo);
                }
        }

	sprintf(gfile->fname,"SL_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_mh%0.1f_rmh%d_MT%d",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->mh,spar->rmigh,MIGRATION_TYPE);
        dnl=200;
	namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;
	char *name=(char *)calloc(namelen,sizeof(char));

        id = (unsigned long)time(NULL);
        
        while(ok==0){
                sprintf(name,"%savInvXplr_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
               
	sprintf(name,"%savInvXplr_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
	
	/****Allocate memory**********************************************/

	avinv=malloc(sizeof(DynVec));
	if (!avinv) { perror("malloc"); exit(1);}
	avinv->sizef=stime->timewindow/stime->tinterval+1;
	avinv->usizef=0;
	avinv->vecf=(double *)calloc(avinv->sizef,sizeof(double));

	/****Dynamics******************/

	spar->plr=0.;
	plrmax=0.1;
	dp=0.01;


	while(spar->plr<=plrmax){
		setCI();

		stime->saveT=stime->transtime;
		stime->Tf=stime->saveT+twind;
		steady=0;
		do{
			callSysDynamics(event,mevent);

			calcRMSError(avinv->vecf,0,avinv->usizef,stats);//test if system reached equilibrium
			if(stats[1]<eps){
				steady=1;
			}else{
				stime->Tf+=twind;
				stime->saveT=stime->Tnow;
				avinv->usizef=0.;
			}
			printf("steady=%d time=%f <avinv>=%f std=%f\n",steady,stime->Tnow,stats[0],stats[1]);
		}while(steady==0);

		nh=listh->usize;
		tot_micr=0.;
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			tot_micr+=spar->micr[idh];
		}
		
		fprintf(gfile->file,"%f %f %f %d %f %f\n",spar->plr,stats[0],stats[1],nh,stime->Tf,tot_micr);
		fflush(gfile->file);
		printf("%f %f %f %d %f %f\n",spar->plr,stats[0],stats[1],nh,stime->Tf,tot_micr);
		
		avinv->usizef=0;

		spar->plr+=dp;
	}
		
	/***freeing memory*****/
	free(avinv->vecf);
	free(avinv);

	if(gfile->file){
		fclose(gfile->file);
		gfile->file=NULL;
	}
	if(gfile->fname){
		free(gfile->fname);
		gfile->fname=NULL;
	}
	if(gfile->fdatapath){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile){
		free(gfile);
		gfile=NULL;
	}
	free(name);
	return;
}
/************************************************************************
*  Store in @SAMPLE files the average investment			*
*  in the system as a function of the host migration coefficient mh     *	
*************************************************************************/
void averInvXmh(Event *event,Event *mevent){
	int i,idh,nh,steady;
        int ok=0,namelen,dnl;
	double eps=0.05,tot_micr,stats[2],mhmax,dmh,aver_rhoe;
	double twind=stime->timewindow;	
	int npar=3,mul,nf;
	double param[npar],expo;
	unsigned long id;
        char nparam[npar][10];
	
       
	/******seting file*************************************/ 
        //generic file struct
	gfile=malloc(sizeof(GenFile));
	if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=300;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");
	
	//file name components
        param[0]=spar->mu;
        param[1]=spar->mig;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
			expo=floor(log10(param[i]));
        		mul=ceil(param[i]/pow(10.,expo));
			nf=0;
			while(ceil(param[i]/pow(10.,expo))!=(int)(param[i]/pow(10.,expo))){
				--expo;
				mul=ceil(param[i]/pow(10.,expo));
				++nf;
			}
                        sprintf(nparam[i],"%de%d",mul,(int)expo);
                }
        }

	#if (NETWORK==0)
	sprintf(gfile->fname,"CG_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1]);
	#else
	sprintf(gfile->fname,"SL_L%d_Ty%d_Kh%d_cost%0.2f_mu%s_mb%s_rmh%d_plr%0.2f_MT%d",L,TYPES,spar->kh,spar->cost,nparam[0],nparam[1],spar->rmigh,spar->plr,MIGRATION_TYPE);
	#endif
        dnl=200;
	namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;
	char *name=(char *)calloc(namelen,sizeof(char));

        id = (unsigned long)time(NULL);
        
        while(ok==0){
                sprintf(name,"%savInvXmh_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
               
	sprintf(name,"%savInvXmh_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
	
	/****Allocate memory**********************************************/

	avinv=malloc(sizeof(DynVec));
	if (!avinv) { perror("malloc"); exit(1);}
	avinv->sizef=stime->timewindow/stime->tinterval+1;
	avinv->usizef=0;
	avinv->vecf=(double *)calloc(avinv->sizef,sizeof(double));

	/****Dynamics******************/

	spar->mh=0.;
	mhmax=10.;
	dmh=0.5;


	while(spar->mh<=mhmax){
		setCI();

		stime->saveT=stime->transtime;
		stime->Tf=stime->saveT+twind;
		steady=0;
		do{
			callSysDynamics(event,mevent);

			calcRMSError(avinv->vecf,0,avinv->usizef,stats);//test if system reached equilibrium
			if(stats[1]<eps){
				steady=1;
			}else{
				stime->Tf+=twind;
				stime->saveT=stime->Tnow;
				avinv->usizef=0.;
			}
			printf("steady=%d time=%f <avinv>=%f std=%f\n",steady,stime->Tnow,stats[0],stats[1]);
		}while(steady==0);

		nh=listh->usize;
		tot_micr=0.;
		aver_rhoe=0.;
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			tot_micr+=spar->micr[idh];
			aver_rhoe+=rho_e[idh];
		}
		aver_rhoe/=(double)nh;
		
		fprintf(gfile->file,"%f %f %f %f %d %f %f\n",spar->mh,stats[0],stats[1],aver_rhoe,nh,stime->Tf,tot_micr);
		fflush(gfile->file);
		printf("mh=%f <w>=%f std(w)=%f <rhoe>=%f nh=%d tf=%f ntot=%f\n",spar->mh,stats[0],stats[1],aver_rhoe,nh,stime->Tf,tot_micr);
		
		avinv->usizef=0;

		spar->mh+=dmh;
	}
		
	/***freeing memory*****/
	free(avinv->vecf);
	free(avinv);

	if(gfile->file){
		fclose(gfile->file);
		gfile->file=NULL;
	}
	if(gfile->fname){
		free(gfile->fname);
		gfile->fname=NULL;
	}
	if(gfile->fdatapath){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile){
		free(gfile);
		gfile=NULL;
	}
	return;
}
/************************************************************************
*               Store heatmap rhXmhX<w>                                 *
*       (cost X bac. migr. rate X average investment)   		*
*************************************************************************/
void costXmbXw(Event *event,Event *mevent){
	int nh,steady;
        int ok=0,namelen,dnl;
	double dmb,dc,cmax,cmin,mbmax,mbmin,eps=0.05,stats[2];
	double twind=stime->timewindow;	
	int mul,nf;
	double expo,param;
	unsigned long id;
        char nparam[10];
	
	/******seting file*************************************/ 
        //generic file struct
	gfile=malloc(sizeof(GenFile));
	if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=300;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");
	
	//file name components
        param=spar->mu;
	
	if(param==0.){
		sprintf(nparam,"0");
	}else{
		expo=floor(log10(param));
		mul=ceil(param/pow(10.,expo));
		nf=0;
		while(ceil(param/pow(10.,expo))!=(int)(param/pow(10.,expo))){
			--expo;
			mul=ceil(param/pow(10.,expo));
			++nf;
		}
		sprintf(nparam,"%de%d",mul,(int)expo);
	}

	#if (NETWORK==0)
	sprintf(gfile->fname,"CG_L%d_Ty%d_Kh%d_mu%s",L,TYPES,spar->kh,nparam);
	#else
	sprintf(gfile->fname,"SL_L%d_Ty%d_Kh%d_mu%s_mh%0.1f_rmh%d_plr%0.2f_MT%d",L,TYPES,spar->kh,nparam,spar->mh,spar->rmigh,spar->plr,MIGRATION_TYPE);
	#endif
        dnl=200;
	namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;
	char *name=(char *)calloc(namelen,sizeof(char));

        id = (unsigned long)time(NULL);
        
        while(ok==0){
		sprintf(name,"%scostXmbXaverW_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
               
	sprintf(name,"%scostXmbXaverW_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
	
	/****Allocate memory**********************************************/

	avinv=malloc(sizeof(DynVec));
	if (!avinv) { perror("malloc"); exit(1);}
	avinv->sizef=stime->timewindow/stime->tinterval+1;
	avinv->usizef=0;
	avinv->vecf=(double *)calloc(avinv->sizef,sizeof(double));

	cmax=0.2;
	cmin=1e-2;
	spar->cost=cmin;
	dc=cmin;
	mbmin=1e-6;
	mbmax=1e-2;

	while(spar->cost<=cmax){
	
		spar->mig=mbmin;
		while(spar->mig<=mbmax){
			setCI();
			stime->saveT=stime->transtime;
			stime->Tf=stime->saveT+twind;
			steady=0;
			do{
				callSysDynamics(event,mevent);
				calcRMSError(avinv->vecf,0,avinv->usizef,stats);//test if system reached equilibrium
				if(stats[1]<eps){
					steady=1;
				}else{
					stime->Tf+=twind;
					stime->saveT=stime->Tnow;
					avinv->usizef=0.;
				}
			
				printf("steady=%d time=%f cost=%f mb=%f <avinv>=%f std=%f\n",steady,stime->Tnow,spar->cost,spar->mig,stats[0],stats[1]);
			}while(steady==0);

			nh=listh->usize;
		
			fprintf(gfile->file,"%f %f %f %f %d %f\n",spar->cost,spar->mig,stats[0],stats[1],nh,stime->Tf);
			fflush(gfile->file);
			printf("%f %f %f %f %d %f\n",spar->cost,spar->mig,stats[0],stats[1],nh,stime->Tf);
		
			avinv->usizef=0;

			if(spar->mig<1e-5){
				dmb=1e-6;
			}else if(spar->mig<1e-4){
				dmb=1e-5;
			}else if(spar->mig<1e-3){
				dmb=1e-4;
			}else if(spar->mig<1e-2){
				dmb=1e-3;
			}else{
				dmb=1e-2;
			}
			spar->mig+=dmb;
		}
		fprintf(gfile->file,"\n");
		printf("\n");

		spar->cost+=dc;
	}
	
	if(gfile->file){
		fclose(gfile->file);
		gfile->file=NULL;
	}
	if(gfile->fname){
		free(gfile->fname);
		gfile->fname=NULL;
	}
	if(gfile->fdatapath){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile){
		free(gfile);
		gfile=NULL;
	}
	
	return;
}
/************************************************************************
*               Store heatmap rhXmhX<w>                                 *
*       (Kh/N X host migration rate coefficient X average investment)   *
*************************************************************************/
void rhXmhXw(Event *event,Event *mevent){
	int i,nh,steady;
        int ok=0,namelen,dnl;
	double dmh,drh,rh,rhmax,mhmax,eps=0.05,stats[2];
	double twind=stime->timewindow;	
	int npar=2,mul,nf;
	double expo,param[npar];
	unsigned long id;
        char nparam[npar][10];
	
       
	/******seting file*************************************/ 
        //generic file struct
	gfile=malloc(sizeof(GenFile));
	if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=300;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");
	
	//file name components
        param[0]=spar->mu;
        param[1]=spar->mig;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
			expo=floor(log10(param[i]));
        		mul=ceil(param[i]/pow(10.,expo));
			nf=0;
			while(ceil(param[i]/pow(10.,expo))!=(int)(param[i]/pow(10.,expo))){
				--expo;
				mul=ceil(param[i]/pow(10.,expo));
				++nf;
			}
                        sprintf(nparam[i],"%de%d",mul,(int)expo);
                }
        }

	#if (NETWORK==0)
	sprintf(gfile->fname,"CG_L%d_Ty%d_cost%0.2f_mu%s_mb%s",L,TYPES,spar->cost,nparam[0],nparam[1]);
	#else
	sprintf(gfile->fname,"SL_L%d_Ty%d_cost%0.2f_mu%s_mb%s_rmh%d_plr%0.2f_MT%d",L,TYPES,spar->cost,nparam[0],nparam[1],spar->rmigh,spar->plr,MIGRATION_TYPE);
	#endif
        dnl=200;
	namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;
	char *name=(char *)calloc(namelen,sizeof(char));

        id = (unsigned long)time(NULL);
        
        while(ok==0){
                sprintf(name,"%srhXmhXaverW_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
               
	sprintf(name,"%srhXmhXaverW_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
	
	/****Allocate memory**********************************************/

	avinv=malloc(sizeof(DynVec));
	if (!avinv) { perror("malloc"); exit(1);}
	avinv->sizef=stime->timewindow/stime->tinterval+1;
	avinv->usizef=0;
	avinv->vecf=(double *)calloc(avinv->sizef,sizeof(double));

	/****Dynamics******************/

	drh=0.05;
	rh=drh;
	spar->kh=N*rh;
	rhmax=0.75;

	mhmax=11.;
	dmh=1.;

	while(rh<=rhmax){
		spar->mh=1.;
		while(spar->mh<mhmax){
			setCI();
			stime->saveT=stime->transtime;
			stime->Tf=stime->saveT+twind;
			steady=0;
		
			do{
				callSysDynamics(event,mevent);
				calcRMSError(avinv->vecf,0,avinv->usizef,stats);//test if system reached equilibrium
				if(stats[1]<eps){
					steady=1;
				}else{
					stime->Tf+=twind;
					stime->saveT=stime->Tnow;
					avinv->usizef=0.;
				}
				printf("steady=%d time=%f <avinv>=%f std=%f\n",steady,stime->Tnow,stats[0],stats[1]);
			}while(steady==0);
			
			nh=listh->usize;
			fprintf(gfile->file,"%f %f %f %f %d %f\n",rh,spar->mh,stats[0],stats[1],nh,stime->Tf);
			fflush(gfile->file);
			printf("rh=%f mh=%f <w>=%f std=%f nh=%d TF=%f\n",rh,spar->mh,stats[0],stats[1],nh,stime->Tf);
			avinv->usizef=0;

			spar->mh+=dmh;
		}
			
		fprintf(gfile->file,"\n");

		rh+=drh;
		spar->kh=N*rh;
	}
		
	/***freeing memory*****/
	free(avinv->vecf);
	free(avinv);

	if(gfile->file){
		fclose(gfile->file);
		gfile->file=NULL;
	}
	if(gfile->fname){
		free(gfile->fname);
		gfile->fname=NULL;
	}
	if(gfile->fdatapath){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile){
		free(gfile);
		gfile=NULL;
	}
	free(name);
	return;
}
/////////////////////////////////////////////////////////////////////////////////////////
/****************************************************************************************
 *		General measuring routines that decide which one to call		*
 ****************************************************************************************/

/***************************************************
*  call routines that measure and store measures   *
*  during the time loop                            *
****************************************************/
void timeMeasures(void){

        #ifdef AVERINVxT
	if((stime->Tnow>=stime->saveT-EPS)&&(stime->Tnow<=stime->saveT+EPS)){
		averInvestmentXt();
		stime->saveT=stime->Tnow+10.;
	}
        #endif
	#ifdef MULTIPLE_COSTS_WxT
	if((stime->Tnow>=stime->saveT-EPS)&&(stime->Tnow<=stime->saveT+EPS)){
		averInvestmentXt();
		stime->saveT=stime->Tnow+10.;
	}
        #endif

        #ifdef SAVE_CONFIG
	if((stime->Tnow>=stime->saveT-EPS)&&(stime->Tnow<=stime->saveT+EPS)){
		printf("Time of measure:%f,  ",stime->Tnow);
		save_config();
		double interval=stime->Tf/(double)NF;
		stime->saveT=stime->Tnow+interval;
		printf("Next time:%f\n",stime->saveT);
	}
        #endif
        #ifdef INV_DIST
	if((stime->Tnow>=stime->saveT-EPS)&&(stime->Tnow<=stime->saveT+EPS)){
		printf("Time of measure:%f,  ",stime->Tnow);
		invDistXt();
		double interval=stime->timewindow/(double)NF;
		stime->saveT=stime->Tnow+interval;
		printf("Next time:%f\n",stime->saveT);
	}
        #endif
	#ifdef NUMHEVENTSxT
	if((stime->Tnow>=stime->saveT-EPS)&&(stime->Tnow<=stime->saveT+EPS)){
		numHostEventsPerDtXt();
		stime->saveT=stime->Tnow+stime->tinterval;
	}
        #endif
	#ifdef CORRxT
	if((stime->Tnow>=stime->saveT-EPS)&&(stime->Tnow<=stime->saveT+EPS)){
		spatialCorrXt();
		stime->saveT=stime->Tnow+stime->tinterval;
	}
        #endif
	#ifdef CLUSTERSxT
	if((stime->Tnow>=stime->saveT-EPS)&&(stime->Tnow<=stime->saveT+EPS)){
		clustersXt();
		stime->saveT=stime->Tnow+stime->tinterval;
	}
        #endif
	#ifdef CLUSTERS_DISTxT
	if((stime->Tnow>=stime->saveT-EPS)&&(stime->Tnow<=stime->saveT+EPS)){
		clustersInvDistXt();
		stime->saveT=stime->Tnow+stime->tinterval;
	}
        #endif
	#ifdef BESTCLUSTER_TIMES
	bestClusterStabilityMeasures();
        #endif
	#ifdef DIFBACOMPxT
	difMicrCompXt();
        #endif

        return;
}

