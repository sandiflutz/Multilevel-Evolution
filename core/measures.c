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

        
	#ifdef GENTIME
	meas->timegh=0.;
	meas->ngh=0;
	meas->nR =(int *)calloc(N,sizeof(int));
	meas->timeR =(double **)calloc(N,sizeof(double *));
	for(i=0; i<N; ++i){
		meas->timeR[i]=(double *)calloc(10,sizeof(int));
	}
	for(i=0; i<N; ++i){
		nR[i]=0;
		memset(meas->timeR[i][nR[i]])=0.;
	}
	#endif
	#ifdef DIFBACOMPxT
	offcomp = malloc(sizeof(DynVec));
        if (!offcomp) { perror("malloc"); exit(1);}
        offcomp->vecf=(double *)calloc(spar->kh,sizeof(double));//max. number of neighbor + focus host
	memset(offcomp->vecf,0.,sizeof(double)*spar->kh);
        offcomp->sizef=spar->kh;
        offcomp->usizef=0;
	#endif

	/**generic file struct to use for storing data during the system evolution**/
	gfile = malloc(sizeof(GenFile));
        if (!gfile) { perror("malloc"); exit(1);}
        gfile->fnsize=400;
        gfile->fname=(char *)calloc(gfile->fnsize,sizeof(char));
	gfile->fdatapath=(char *)malloc(sizeof(char)*50);
        sprintf(gfile->fdatapath,"data_manipulation/");

	/**file name components**/

        char *ngeral=(char *)calloc(250,sizeof(char));
        char *ntneg = (char *)calloc(50,sizeof(char));
	int npar=5;
        char nparam[npar][10];
        double param[npar];
        param[0]=Bacv;
        param[1]=spar->cost;
        param[2]=spar->mu;
        param[3]=spar->mig;
        param[4]=Fmin;
        
	for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
                        sprintf(nparam[i],"1e%d",(int)log10(param[i]));
                }
        }
	#if(TV==1)
        sprintf(ngeral,"N%d_Ty%d_Kh%d_net%d_Gh%d_CI%d_TV%d_Bv%s_Fmin%s_cost%s_mu%s_mb%s_mh%0.1f",N,TYPES,spar->kh,NETWORK,Gh,CI,TV,nparam[0],nparam[4],nparam[1],nparam[2],nparam[3],spar->mh);
	#else
        sprintf(ngeral,"N%d_Ty%d_Kh%d_net%d_Gh%d_CI%d_TV%d_Bv%s_cost%s_mu%s_mb%s_mh%0.1f",N,TYPES,spar->kh,NETWORK,Gh,CI,TV,nparam[0],nparam[1],nparam[2],nparam[3],spar->mh);
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
	
	if(gfile->file!=NULL){
		fclose(gfile->file);
		gfile->file = NULL;
	}
	if(gfile->fname!=NULL){
		free(gfile->fname);
		gfile->fname = NULL;
	}
	if(gfile->fdatapath!=NULL){
		free(gfile->fdatapath);
		gfile->fdatapath=NULL;
	}
	if(gfile!=NULL){
		free(gfile);
		gfile=NULL;
	}
	#ifdef GENTIME
	for(i=0; i<N; ++i){
		free(timeR[i]);
	}
	free(timeR);
	free(nR);
	#endif
	#ifdef DIFBACOMPxT
	if(offcomp->vecf!=NULL)free(offcomp->vecf);
        if(offcomp!=NULL)free(offcomp);
	#endif

	return;
}
/******************************************
*	Open Global Files                 *
*******************************************/
void openFiles(void){
#ifdef TMEAS
        int ok=0,namelen,dnl;
        unsigned long id;

        id = (unsigned long)time(NULL);

        dnl=200;
        namelen=strlen(gfile->fname)+strlen(gfile->fdatapath)+dnl;

        char *name=(char *)calloc(namelen,sizeof(char));

#ifdef DENSb1xT
        /*each execution will produce a file with a different name (with a "random" id at the end of the name)*/
        while(ok==0){
                sprintf(name,"%sdensb1Xt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
        sprintf(name,"%sdensb1Xt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}

#endif
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
#ifdef EmptyFreqxT
        while(ok==0){
                sprintf(name,"%semptyFreqXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
	sprintf(name,"%semptyFreqXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
        gfile->file=fopen(name,"w");
        if (gfile->file==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef MEANBFRACxT
        while(ok==0){
                sprintf(name,"%smeanFracBXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
        sprintf(name,"%smeanFracBXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
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
#ifdef GENTIME
        while(ok==0){
                sprintf(name,"%sgentimeXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
                gfile->file=fopen(name,"r");
                if(gfile->file!=NULL){
                        ++id;
                        fclose(gfile->file);
                }else{
                        ok=1;
                }
        }
	sprintf(name,"%sgentimeXt_%s_%ld.dat",gfile->fdatapath,gfile->fname,id);
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
			densInvH[i]+=bac[idh][j]*spar->inv[j];
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
			freqInvH[i]+=bac[idh][j]*spar->inv[j]/spar->micr[idh];
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
/********************************************************
*  Calculate the difference in microbial               *
*  composition between host @idp and its children @idk.*
*  Store result in a vector where each element         *
*  correspond to a different event                     *
********************************************************/
void storeBacDiffComp(int idp,int idk){
	int j,new;
        double mdiff;

        mdiff=0.;
        new=0;
        for(j=0; j<TYPES; ++j){
                if(bac[idp][j]!=0.){
			mdiff+=spar->inv[j]*(bac[idp][j]/spar->micr[idp]-bac[idk][j]/spar->micr[idk]);
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
/************************************************
*   store the current average investment level  *
*   in an isolate host                          *
**************************************************/
void densB1Xt(void){
	int j;
	double averinv;

	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:time | 2:average investment level on an isolated host | 5:abundance of the micr. pop. | 6:#of time steps)\n");
	}
	
		
	averinv=0.;
	for(j=0; j<TYPES; ++j){
		averinv+=bac[0][j]*spar->inv[j];
	}
	averinv/=spar->micr[0];
	
	fprintf(gfile->file,"%f %f %f\n",stime->Tnow,averinv,spar->micr[0]);
	printf("time=%f averinv=%f micr=%f\n",stime->Tnow,averinv,spar->micr[0]);
        
	return;
}
/********************************************************
* Stores the mean fraction of each type of bacteria in  *
* the system                                            *
*********************************************************/
void meanFracXt(void){
	int i,j,idh,nh;
	double *avbacfreq;
	
	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:time 2:mean freq. of bac of type=%d ",TYPES-1);
		for(j=TYPES-2; j>=0; --j){
			fprintf(gfile->file,"%d:mean freq. of bac of type=%d ",TYPES-j+1,j);
		}
		fprintf(gfile->file,"%d:numsteps %d:nh\n",TYPES+2,TYPES+3);
	}
	
	nh=listh->usize;

	avbacfreq=(double *)calloc(TYPES,sizeof(double));
	memset(avbacfreq,0.,sizeof(double));
	
        fprintf(gfile->file,"%f ",stime->Tnow);
        printf("%f ",stime->Tnow);
	for(j=TYPES-1; j>=0; --j){
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			avbacfreq[j]+=bac[idh][j]/(spar->micr[idh]);
		}
		avbacfreq[j]/=(double)nh;
		/*storing data*/
		fprintf(gfile->file,"%f ",avbacfreq[j]);
		printf("<bac[%d]>=%f ",j,avbacfreq[j]);
	}
        fprintf(gfile->file,"%d\n",nh);
        printf("nh=%d\n",nh);

	free(avbacfreq);
	return;
}
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
	int i,idh,nh;
	double averinv,tot_micr;
	
	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:time 2:average cumulative investment 3:average microbial density 4:#of hosts 5:dens_host 6:#of time steps 7:dth\n");
	}

	/****setting investiment density per host vector and calculating average investment in the system*****/
	nh=listh->usize;
	averinv=calcAverInv();

	tot_micr=0.;
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		tot_micr+=spar->micr[idh];
	}

	/*storing data*/
        fprintf(gfile->file,"%f %f %f %d %f %f\n",stime->Tnow,averinv,(double)tot_micr/nh,nh,(double)nh/N,stime->dth);
        printf("t=%f averinv=%f avmicrdens=%f nh=%d nh/N=%f dth=%f\n",stime->Tnow,averinv,(double)tot_micr/nh,nh,(double)nh/N,stime->dth);

        return;
}
/********************************************************
*	Store measures related to vancancy		*
*	frequency as functions	of time			*
*	-> average group vacancy freq.,<rho_e>		*
*	-> standart deviation of rho_e			*
*	-> (Kh-H(t))/N=Rho_e-Rho_ekh			*
*	(where Rho_e is the system vacancy 		*
*	freq. and Rho_ekh is the vacancy freq.		*
*	related to the carrying capacity=(N-Kh)/N	*
*********************************************************/
void emptyFreqXt(void){
	int i,idh,nh;
	double averinv,stats[2],*eff_rhoe;

	nh=listh->usize;
	eff_rhoe=(double *)calloc(nh,sizeof(double));
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		eff_rhoe[i]=rho_e[idh];
	}
	
	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:time 2:<rho_e> 3:std(rho_e) 4:Rho_e-Rho_ekh 5:<w>\n");
	}
	averinv=calcAverInv();

	calcRMSError(eff_rhoe,0,nh-1,stats);
	
	/*storing data*/
        fprintf(gfile->file,"%f %f %f %f %f\n",stime->Tnow,stats[0],stats[1],averinv,(double)(spar->kh-nh)/N);
        printf("t=%f <rho_e>=%f std=%f <w>=%f (Kh-nh)/N=%f\n",stime->Tnow,stats[0],stats[1],averinv,(double)(spar->kh-nh)/N);

	free(eff_rhoe);
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
                pointsize=0.5;
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
        sprintf(nt_format,"Tf%d",stime->Tnow);
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
        fprintf(fgp,"set output'%s_Tf%d_fm.eps'\n",name,stime->Tnow);
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
/***************************************************
*  calculate spatial corretation: both horizontal  *
*  and vertical for a specific distance,           *
*  considering                                     *
****************************************************/
void spatialCorrXt(void){
	int dist;
	double corr_tot,corr[2];

	if(stime->Tnow==0.){
		if(gfile->file==NULL){
			printf("You are trying to write in a file that doesn't exist.\n");
			exit(1);
		}else{
			fprintf(gfile->file,"#1:time 2:corrx 3:corry 4:corr 5:#of hosts\n");
			printf("#1:time 2:corrx 3:corry 4:corr 5:#of hosts\n");
		}
	}

	dist=1;
	corr_tot=spatialCorr(host,N,dist,RIGHT,DOWN,neighbor,corr);/*sending: 1-state vector,2-square lattice size, 3-distance for calculating spatial correlation
						 *4-index of horizontal neighbors (right or left), 5-index of vertical neighbors (top or bottom)
						 *5-vector for storing vertical and horizontal correlations*/ 
	fprintf(gfile->file,"%f %f %f %f %d\n",stime->Tnow,corr[0],corr[1],corr_tot,listh->usize);
	printf("t=%f corrx=%f corry=%f corr=%f nh=%d\n",stime->Tnow,corr[0],corr[1],corr_tot,listh->usize);

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
/***************************************************
*  storing average host generation time in a file  *
****************************************************/
void genHostTime(void){
	double mean_gtime;
	if(stime->Tnow==0.){
		fprintf(gfile->file,"#1:time 2:#mean host generation time 3:#of repr. events used 4:#of hosts\n");
	}
	
	if(meas->ngh==0){
		mean_gtime=0.;
	}else{
		mean_gtime=(double)meas->timegh/meas->ngh;
	}
	fprintf(gfile->file,"%f %f %d %d\n",stime->Tnow,mean_gtime,meas->ngh,listh->usize);
	printf("t=%f <timegh>=%f sample=%d nh=%d\n",stime->Tnow,mean_gtime,meas->ngh,listh->usize);

	return;
}
/*******steady state measures***********************************************************/

/****************************************************************
*  Store in @SAMPLE files the average investment		*
*  in the system as a function of the system carrying dilution	*	
*****************************************************************/
void averInvXrh(Event *event,Event *mevent){
	int i,idh,nh,steady;
        int ok=0,namelen,dnl;
	double rh,drh=0.05,eps=0.05,tot_micr,stats[2];
	double twind=stime->timewindow;	
	int npar=5;
	double param[npar];
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
        param[0]=Bacv;
        param[1]=spar->cost;
        param[2]=spar->mu;
        param[3]=spar->mig;
	param[4]=Fmin;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
                        sprintf(nparam[i],"1e%d",(int)log10(param[i]));
                }
        }

	#if (TV==1)
	sprintf(gfile->fname,"N%d_Ty%d_net%d_Gh%d_CI%d_TV%d_Bv%s_Fmin%s_cost%s_mu%s_mb%s_mh%0.1f",N,TYPES,NETWORK,Gh,CI,TV,nparam[0],nparam[4],nparam[1],nparam[2],nparam[3],spar->mh);
	#else
	sprintf(gfile->fname,"N%d_Ty%d_net%d_Gh%d_CI%d_TV%d_Bv%s_cost%s_mu%s_mb%s_mh%0.1f",N,TYPES,NETWORK,Gh,CI,TV,nparam[0],nparam[1],nparam[2],nparam[3],spar->mh);
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

	fclose(gfile->file);
	free(gfile->fname);
	free(gfile->fdatapath);
	free(gfile);
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
	int npar=5;
	double param[npar];
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
        param[0]=Bacv;
        param[1]=spar->cost;
        param[2]=spar->mu;
        param[3]=spar->mig;
	param[4]=Fmin;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
                        sprintf(nparam[i],"1e%d",(int)log10(param[i]));
                }
        }

	#if (TV==1)
	sprintf(gfile->fname,"N%d_Ty%d_net%d_Kh%d_CI%d_TV%d_Bv%s_Fmin%s_cost%s_mu%s_mb%s_mh%0.1f",N,TYPES,NETWORK,spar->kh,CI,TV,nparam[0],nparam[4],nparam[1],nparam[2],nparam[3],spar->mh);
	#else
	sprintf(gfile->fname,"N%d_Ty%d_net%d_Kh%d_CI%d_TV%d_Bv%s_cost%s_mu%s_mb%s_mh%0.1f",N,TYPES,NETWORK,spar->kh,CI,TV,nparam[0],nparam[1],nparam[2],nparam[3],spar->mh);
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

	fclose(gfile->file);
	free(gfile->fname);
	free(gfile->fdatapath);
	free(gfile);
	free(name);
	return;
}
/************************************************************************
*  Store in @SAMPLE files the average investment			*
*  in the system as a function of the microbial migration rate		*	
*************************************************************************/
void averInvXmb(Event *event,Event *mevent){
	int i,idh,nh,steady;
        int ok=0,namelen,dnl;
	double dmb,mbmax,eps=0.05,tot_micr,stats[2];
	double twind=stime->timewindow;	
	int npar=5;
	double param[npar];
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
        param[0]=Bacv;
        param[1]=spar->cost;
        param[2]=spar->mu;
	param[3]=Fmin;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
                        sprintf(nparam[i],"1e%d",(int)log10(param[i]));
                }
        }

	#if (TV==1)
	sprintf(gfile->fname,"N%d_Ty%d_net%d_Kh%d_Gh%d_CI%d_TV%d_Bv%s_Fmin%s_cost%s_mu%s_mh%0.1f",N,TYPES,NETWORK,spar->kh,spar->gh,CI,TV,nparam[0],nparam[3],nparam[1],nparam[2],spar->mh);
	#else
	sprintf(gfile->fname,"N%d_Ty%d_net%d_Kh%d_Gh%d_CI%d_TV%d_Bv%s_cost%s_mu%s_mh%0.1f",N,TYPES,NETWORK,spar->kh,spar->gh,CI,TV,nparam[0],nparam[1],nparam[2],spar->mh);
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

	spar->mig=1e-6;
	mbmax=0.1;
	dmb=1e-3;

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

		spar->mig+=dmb;
	}
		
	/***freeing memory*****/
	free(avinv->vecf);
	free(avinv);

	fclose(gfile->file);
	free(gfile->fname);
	free(gfile->fdatapath);
	free(gfile);
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
	int npar=5;
	double param[npar];
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
        param[0]=Bacv;
        param[1]=spar->mu;
        param[2]=spar->mig;
	param[3]=Fmin;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
                        sprintf(nparam[i],"1e%d",(int)log10(param[i]));
                }
        }

	#if (TV==1)
	sprintf(gfile->fname,"N%d_Ty%d_net%d_Kh%d_Gh%d_CI%d_TV%d_Bv%s_Fmin%s_mu%s_mb%s_mh%0.1f",N,TYPES,NETWORK,spar->kh,spar->gh,CI,TV,nparam[0],nparam[3],nparam[1],nparam[2],spar->mh);
	#else
	sprintf(gfile->fname,"N%d_Ty%d_net%d_Kh%d_Gh%d_CI%d_TV%d_Bv%s_mu%s_mb%s_mh%0.1f",N,TYPES,NETWORK,spar->kh,spar->gh,CI,TV,nparam[0],nparam[1],nparam[2],spar->mh);
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
	cost_max=1.;
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

	fclose(gfile->file);
	free(gfile->fname);
	free(gfile->fdatapath);
	free(gfile);
	free(name);
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
	int npar=5;
	double param[npar];
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
        param[0]=Bacv;
	param[1]=spar->cost;
        param[2]=spar->mu;
        param[3]=spar->mig;
	param[4]=Fmin;

        for(i=0; i<npar; ++i){
                if(param[i]==0.){
                        sprintf(nparam[i],"0");
                }else{
                        sprintf(nparam[i],"1e%d",(int)log10(param[i]));
                }
        }

	#if (TV==1)
	sprintf(gfile->fname,"N%d_Ty%d_net%d_Gh%d_CI%d_TV%d_Bv%s_Fmin%s_cost%s_mu%s_mb%s",N,TYPES,NETWORK,spar->gh,CI,TV,nparam[0],nparam[4],nparam[1],nparam[2],nparam[3]);
	#else
	sprintf(gfile->fname,"N%d_Ty%d_net%d_Gh%d_CI%d_TV%d_Bv%s_cost%s_mu%s_mb%s",N,TYPES,NETWORK,spar->gh,CI,TV,nparam[0],nparam[1],nparam[2],nparam[3]);
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

	fclose(gfile->file);
	free(gfile->fname);
	free(gfile->fdatapath);
	free(gfile);
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

        #ifdef DENSb1xT
	if((stime->Tnow<=stime->Tf)){
		densB1Xt();
	}
        #endif
        #ifdef AVERINVxT
//	if((stime->Tnow>=stime->saveT-err)&&(stime->Tnow<=stime->saveT+err)){
		averInvestmentXt();
//		stime->saveT+=stime->tinterval;
//	}
        #endif
        #ifdef EmptyFreqxT
		 emptyFreqXt();
        #endif
        #ifdef MEANBFRACxT
	if((stime->Tnow<=stime->Tf)){
		meanFracXt();
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
		stime->saveT+=stime->tinterval;
	}
        #endif
	#ifdef CORRxT
	if((stime->Tnow>=stime->saveT-EPS)&&(stime->Tnow<=stime->saveT+EPS)){
		spatialCorrXt();
		stime->saveT+=stime->tinterval;
	}
        #endif
	#ifdef GENTIME
	genHostTime();
        #endif
	#ifdef DIFBACOMPxT
	difMicrCompXt();
        #endif

        return;
}

