/* meaures.c */
#include <time.h>
#include"globals.h"
#include"tools.h"
#include"measures.h"
/***********************************************************
*     allocate memory for global arrays and structs        *
*     related to time measuraments                         *
************************************************************/
void allocateMemTM(TimeMeasures *meas){
	int i;
        char *ngeral=(char *)calloc(200,sizeof(char));
        char *nevo = (char *)calloc(50,sizeof(char));
        char *ntneg = (char *)calloc(50,sizeof(char));
        char nparam[4][10];
        double param[4];

	meas->idh_h1=0.;
        meas->Ti=0.;
        meas->Tf=TF;
        meas->NTf=NTf_me;
	meas->dth=Dt_ref;
        
	meas->saveT=0.;
        meas->nfiles=0;
        
	meas->timegh=0.;
	meas->ngh=0;
	#ifdef GENTIME
	timeb =(double *)calloc(N,sizeof(double));
	for(i=0; i<N; ++i){
		if(host[i]==1){
			timeb[i]=0.;
		}else{
			timeb[i]=-1.;
		}
	}
	#endif

        meas->ftnpars_size=400;
        meas->ftname_pars=(char *)calloc(meas->ftnpars_size,sizeof(char));

        param[0]=Bacv;
        param[1]=spar->cost;
        param[2]=spar->mu;
        param[3]=spar->mig;
        
	for(i=0; i<4; ++i){
                if(param[i]==0){
                        sprintf(nparam[i],"0");
                }else{
                        sprintf(nparam[i],"1e%d",(int)log10(param[i]));
                }
        }
        sprintf(ngeral,"N%d_Ty%d_Tp%d_Tn%d_Kh%d_net%d_Gh%d_CI%d_sigma%0.2f_Bv%s_cost%s_mu%s_mig%s",N,TYPES,Tpos,Tneg,spar->kh,NETWORK,Gh,CI,spar->sigma,nparam[0],nparam[1],nparam[2],nparam[3]);
	#if (EVO==0)
	sprintf(nevo,"_v0");
	#elif (EVO==1)
	sprintf(nevo,"_tlp");
	#else
	sprintf(nevo,"_mcs");
	#endif
                
	#if (Tneg>0)
	sprintf(ntneg,"_CRnnA%d_CRnnB%0.1f_CRnpA%d_CRnpB%0.1f",(int)CRnn0,CRnn1,(int)CRnp0,CRnp1);
	#else
	sprintf(ntneg,"");
	#endif
        sprintf(meas->ftname_pars,"%s%s%s",ngeral,nevo,ntneg);
        free(ngeral);
        free(nevo);
        free(ntneg);

	return;
}
/************************************************
*       Free Allocated Memory for time		*
*       measurements				*
*************************************************/
void freeMemTM(TimeMeasures *meas){
	
	free(meas->ftname_pars);
	#ifdef GENTIME
	free(timeb);
	#endif

	return;
}
/******************************************
*	Open Global Files                 *
*******************************************/
void openFiles(TimeMeasures *meas){
        int ok=0,namelen,dnl;
        unsigned long id;

        id = (unsigned long)time(NULL);

        dnl=200;
        namelen=strlen(meas->ftname_pars)+strlen(fdatapath)+dnl;

        char *name=(char *)calloc(namelen,sizeof(char));

#ifdef DENSb1xT
        /*each execution will produce a file with a different name (with a "random" id at the end of the name)*/
        while(ok==0){
                sprintf(name,"%sdensb1Xt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
                meas->file_tmeas=fopen(name,"r");
                if(meas->file_tmeas!=NULL){
                        ++id;
                        fclose(meas->file_tmeas);
                }else{
                        ok=1;
                }
        }
        sprintf(name,"%sdensb1Xt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
        meas->file_tmeas=fopen(name,"w");
        if (meas->file_tmeas==NULL) { perror("malloc"); exit(1);}

#endif
#ifdef AVERINVxT
        while(ok==0){
                sprintf(name,"%saverInvXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
                meas->file_tmeas=fopen(name,"r");
                if(meas->file_tmeas!=NULL){
                        ++id;
                        fclose(meas->file_tmeas);
                }else{
                        ok=1;
                }
        }
        sprintf(name,"%saverInvXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
        meas->file_tmeas=fopen(name,"w");
        if (meas->file_tmeas==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef MEANBFRACxT
        while(ok==0){
                sprintf(name,"%smeanFracBXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
                meas->file_tmeas=fopen(name,"r");
                if(meas->file_tmeas!=NULL){
                        ++id;
                        fclose(meas->file_tmeas);
                }else{
                        ok=1;
                }
        }
        sprintf(name,"%smeanFracBXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
        meas->file_tmeas=fopen(name,"w");
        if (meas->file_tmeas==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef NUMHEVENTSxT
        while(ok==0){
                sprintf(name,"%snumheventsXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
                meas->file_tmeas=fopen(name,"r");
                if(meas->file_tmeas!=NULL){
                        ++id;
                        fclose(meas->file_tmeas);
                }else{
                        ok=1;
                }
        }
        sprintf(name,"%snumheventsXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
        meas->file_tmeas=fopen(name,"w");
        if (meas->file_tmeas==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef CORRxT
        while(ok==0){
                sprintf(name,"%scorrXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
                meas->file_tmeas=fopen(name,"r");
                if(meas->file_tmeas!=NULL){
                        ++id;
                        fclose(meas->file_tmeas);
                }else{
                        ok=1;
                }
        }
	sprintf(name,"%scorrXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
        meas->file_tmeas=fopen(name,"w");
        if (meas->file_tmeas==NULL) { perror("malloc"); exit(1);}
#endif
#ifdef GENTIME
        while(ok==0){
                sprintf(name,"%sgentimeXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
                meas->file_tmeas=fopen(name,"r");
                if(meas->file_tmeas!=NULL){
                        ++id;
                        fclose(meas->file_tmeas);
                }else{
                        ok=1;
                }
        }
	sprintf(name,"%sgentimeXt_%s_%ld.dat",fdatapath,meas->ftname_pars,id);
        meas->file_tmeas=fopen(name,"w");
        if (meas->file_tmeas==NULL) { perror("malloc"); exit(1);}
#endif
        free(name);
        return;
}
/**********************************************
*   		Close Global Files            *                    
***********************************************/
void closeFiles(TimeMeasures *meas){
	
	if(meas->file_tmeas!=NULL){
		fclose(meas->file_tmeas);
	}

	return;
}
/**********************************************
 * calculate investment density per Host     *
 **********************************************/
void calcInvDens(double *densInvH){
	int i,idh,j,nh;
	
	nh=listh->usize;
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		densInvH[i]=0.;
		for(j=0; j<TYPES; ++j){
			densInvH[i]+=bac[idh][j]*spar->s[j]*spar->inv[j];
		}
	}

	return;
}
/**********************************************
 * calculate investment frequency per Host    *
 **********************************************/
void calcInvFreq(double *freqInvH){
	int i,idh,j,nh;
	
	nh=listh->usize;
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		freqInvH[i]=0.;
		for(j=0; j<TYPES; ++j){
			freqInvH[i]+=bac[idh][j]*spar->s[j]*spar->inv[j]/spar->micr[idh];
		}
	}

	return;
}
/*******************************************************************
 * calculate investment distribution among                         *
 * hosts                                                           *
 *******************************************************************/
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
/************************************************
*   store the current average investment level  *
*   in an isolate host                          *
**************************************************/
void densB1Xt(TimeMeasures *meas){
	int j;
	double averinv;

	if(meas->Tnow==meas->Ti){
		fprintf(meas->file_tmeas,"#1:time | 2:average investment level on an isolated host | 3:frac. of bac. of type 0 | 4:frac. of bac. of type 1 | 5:abundance of the micr. pop. | 6:#of time steps)\n");
	}
	
	averinv=0.;
	for(j=0; j<TYPES; ++j){
		averinv+=bac[meas->idh_h1][j]*spar->s[j]*spar->inv[j]/spar->micr[meas->idh_h1];
	}
	
	fprintf(meas->file_tmeas,"%f %f %f %f %f %d\n",meas->Tnow,averinv,bac[meas->idh_h1][0]/spar->micr[meas->idh_h1],bac[meas->idh_h1][1]/spar->micr[meas->idh_h1],spar->micr[meas->idh_h1],meas->NTnow);
	printf("time=%f averinv=%f fbac0=%f fbac1=%f micr=%f numsteps=%d\n",meas->Tnow,averinv,bac[meas->idh_h1][0]/spar->micr[meas->idh_h1],bac[meas->idh_h1][1]/spar->micr[meas->idh_h1],spar->micr[meas->idh_h1],meas->NTnow);
       	printf("%f %f\n",meas->Tnow,averinv);
        
	return;
}
/********************************************************
* Stores the mean fraction of each type of bacteria in  *
* the system                                            *
*********************************************************/
void meanFracXt(TimeMeasures *meas){
	int i,j,idh,nh;
	double *avbacfreq;
	
	if(meas->Tnow==meas->Ti){
		fprintf(meas->file_tmeas,"#1:time 2:mean freq. of bac of type=%d ",TYPES-1);
		for(j=TYPES-2; j>=0; --j){
			fprintf(meas->file_tmeas,"%d:mean freq. of bac of type=%d ",TYPES-j+1,j);
		}
		fprintf(meas->file_tmeas,"%d:numsteps\n",TYPES+2);
	}
	
	nh=listh->usize;

	avbacfreq=(double *)calloc(TYPES,sizeof(double));
	memset(avbacfreq,0.,sizeof(double));
	
        fprintf(meas->file_tmeas,"%f ",meas->Tnow);
        printf("%f ",meas->Tnow);
	for(j=TYPES-1; j>=0; --j){
		for(i=0; i<nh; ++i){
			idh=listh->vec[i];
			avbacfreq[j]+=bac[idh][j]/(spar->micr[idh]);
		}
		avbacfreq[j]/=(double)nh;
		/*storing data*/
		fprintf(meas->file_tmeas,"%f ",avbacfreq[j]);
		printf("%f ",avbacfreq[j]);
	}
        fprintf(meas->file_tmeas,"%d\n",meas->NTnow);
        printf("%d\n",meas->NTnow);

	free(avbacfreq);
	return;
}
/********************************************************
*   stores the average investment in the system         *
*	-Def: avInv=sum_ij(bac_ij*inv_j)/sum_ij(bacij)  *
*	where i is the host index, j is the type of     *
*   	bacteria index, bac_ij is the abundance of      *
*   	type j in host i and inv_j is the investment    *
*   	of type j.                                      *
*   	-When TYPES=2:                                  *
*   		-inv is either 0 (for neutrals)         *
*   		or 1 (for helpers)                      *
*   		-avInv is the mean frequency of         *
*   		helpers in the system                   *
*********************************************************/
void averInvestmentXt(TimeMeasures *meas){
	int i,idh,nh;
	double averinv,*densInvH,tot_micr;
	
	if(meas->Tnow==meas->Ti){
		fprintf(meas->file_tmeas,"#1:time 2:average cumulative investment 3:average microbial density 4:#of hosts 5:dens_host 6:#of time steps\n");
	}

	/****setting investiment density per host vector*****/
	nh=listh->usize;
	densInvH=(double *)calloc(nh,sizeof(double));/*indexes are the same as in @listh->vec[] 
						      * (real host indexes are idh=listh->vec[i] (for 0<=i<nh)) */
	calcInvDens(densInvH);

	/*********calculating average investment is the system***********/
	averinv=0.;
	tot_micr=0.;
        for(i=0; i<nh; ++i){
                idh=listh->vec[i];
		tot_micr+=spar->micr[idh];
		averinv+=densInvH[i];
        }
	averinv/=tot_micr;//averaging over the whole microbe population
	

	/*storing data*/
        fprintf(meas->file_tmeas,"%f %f %f %d %f %d\n",meas->Tnow,averinv,(double)tot_micr/nh,nh,(double)nh/N,meas->NTnow);
        printf("t=%f averinv=%f avmicrdens=%f nh=%d nh/N=%f numsteps=%d\n",meas->Tnow,averinv,(double)tot_micr/nh,nh,(double)nh/N,meas->NTnow);

	/*freeing allocated memory*/
	free(densInvH);
        return;
}
/**************************************
*  snapshot of host network           *
*  The colors indicate the acumulated *
*  investment of each host            *
***************************************/
void save_config(TimeMeasures *meas){
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
	dnl=100;
        cnamelen=strlen(meas->ftname_pars)+dnl;//core name length
        
	//core name structure
        name=(char *)calloc(cnamelen,sizeof(char));//base name
        sprintf(name,"snapshot_%s",meas->ftname_pars);
	
	nt_format=(char *)calloc(dnl,sizeof(char));//format of the time term, on the name (files for animations have to be ordered: better to use a 0<t<1 format)
#if (FIG_EXT==0)
        sprintf(nt_format,"idt%f",meas->Tnow/meas->Tf);
#else
        sprintf(nt_format,"Tf%f",meas->Tnow);
#endif
	
	//file names for data and for gnuplot script
	snamelen=strlen(name)+strlen(fdatapath)+strlen(nt_format)+dnl;//length of the specific file names
	namedat=(char *)calloc(snamelen,sizeof(char));//data file name
        name_gp=(char *)calloc(snamelen,sizeof(char));//name for the gnuplot data file


        sprintf(namedat,"%s%s_%s.dat",fdatapath,name,nt_format);
        sprintf(name_gp,"%s%s_%s.gp",fdatapath,name,nt_format);

        fconfig = fopen(namedat,"w");
        fgp = fopen(name_gp,"w");


	/*********Defining point size for the gnuplot script****************************************/
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
	fprintf(fgp,"set output'%s_idt%f.png'\n",name,meas->Tnow/meas->Tf);
#else
        fprintf(fgp,"set term post eps enha color 20\n");
        fprintf(fgp,"set output'%s_Tf%f_fm.eps'\n",name,meas->Tnow);
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
 * store investment distribution among hosts and     *
 * create a gnuplot script to create graphics        *
 *****************************************************/
void invDistXt(TimeMeasures *meas){
	int i,nh,nbins=100,namelen,dnl;
	double binsize,*hist_inv,*freqInvH;
        char *name,*namedat,*name_gp,*nt_format;
        FILE *fhist,*fgp;

	/****creating files (if there are file with the same names, there are subscribed)********/
	dnl=100;
	namelen=strlen(meas->ftname_pars)+strlen(fdatapath)+dnl;
        nt_format=(char *)calloc(dnl,sizeof(char));//base name
        name=(char *)calloc(namelen,sizeof(char));//format of the time term, on the name (files for animations have to be ordered: better to use a 0<t<1 format)
        namedat=(char *)calloc((namelen+dnl),sizeof(char));//data file name
        name_gp=(char *)calloc((namelen+dnl),sizeof(char));//name for the gnuplot data file
        
	sprintf(name,"invDistXt_%s",meas->ftname_pars);

#if (FIG_EXT==0)
        sprintf(nt_format,"idt%f",(double)meas->Tnow/meas->Tf);
#else
        sprintf(nt_format,"Tf%d",meas->Tnow);
#endif
        sprintf(namedat,"%s%s_%s.dat",fdatapath,name,nt_format);
        sprintf(name_gp,"%s%s_%s.gp",fdatapath,name,nt_format);
        
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
	
	if(meas->Tnow==meas->Ti)fprintf(fhist,"#1:average inv. in hosts 2:frac. of hosts 3:numsteps\n");
	for(i=0; i<nbins; ++i){
		fprintf(fhist,"%f %f %f %f\n",(double)i*binsize,(double)hist_inv[i]/nh,freqInvH[0],meas->Tnow);
	}
	
	/*******************filling gnuplot file*****************************/

#if (FIG_EXT==0)
       fprintf(fgp,"set term png size 720,540\n");
       fprintf(fgp,"set output'%s_idt%f.png'\n",name,(double)meas->Tnow/meas->Tf);
#else
        fprintf(fgp,"set term post eps enha color 20\n");
        fprintf(fgp,"set output'%s_Tf%d_fm.eps'\n",name,meas->Tnow);
#endif

	fprintf(fgp,"set title'{/=15 time steps=%f}'\n",meas->Tnow);
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
void numHostEventsPerDtXt(TimeMeasures *meas){

	if(meas->Tnow==meas->Ti){
		fprintf(meas->file_tmeas,"#1:time 2:#of host events per Dt_ref 3:#of host births per Dt_ref 4:#of host deaths per Dt_ref 5:#of hosts\n");
	}
	fprintf(meas->file_tmeas,"%f %d %d %d %d\n",meas->Tnow,meas->numb+meas->numd,meas->numb,meas->numd,listh->usize);
	printf("%f %d %d %d %d\n",meas->Tnow,meas->numb+meas->numd,meas->numb,meas->numd,listh->usize);

	return;
}
/***************************************************
*  calculate spatial corretation: both horizontal  *
*  and vertical for a specific distance,           *
*  considering                                     *
****************************************************/
void spatialCorrXt(TimeMeasures *meas){
	int dist;
	double corr_tot,corr[2];

	if(meas->Tnow==meas->Ti){
		if(meas->file_tmeas==NULL){
			printf("You are trying to write in a file that doesn't exist.\n");
			exit(1);
		}else{
			fprintf(meas->file_tmeas,"#1:time 2:corrx 3:corry 4:corr 5:#of hosts\n");
			printf("#1:time 2:corrx 3:corry 4:corr 5:#of hosts\n");
		}
	}

	dist=1;
	corr_tot=spatialCorr(host,N,dist,RIGHT,DOWN,neighbor,&corr);/*sending: 1-state vector,2-square lattice size, 3-distance for calculating spatial correlation
						 *4-index of horizontal neighbors (right or left), 5-index of vertical neighbors (top or bottom)
						 *5-vector for storing vertical and horizontal correlations*/ 
	fprintf(meas->file_tmeas,"%f %f %f %f %d\n",meas->Tnow,corr[0],corr[1],corr_tot,listh->usize);
	printf("t=%f corrx=%f corry=%f corr=%f nh=%d\n",meas->Tnow,corr[0],corr[1],corr_tot,listh->usize);

	return;
}
/***************************************************
*  calculating average host generation time        *
****************************************************/
void genHostTime(TimeMeasures *meas){
	double mean_gtime;
	if(meas->Tnow==meas->Ti){
		fprintf(meas->file_tmeas,"#1:time 2:#mean host generation time 3:#of repr. events used 4:#of hosts\n");
	}
	
	if(meas->ngh==0){
		mean_gtime=0.;
	}else{
		mean_gtime=(double)meas->timegh/meas->ngh;
	}
	fprintf(meas->file_tmeas,"%f %f %d %d\n",meas->Tnow,mean_gtime,meas->ngh,listh->usize);
	printf("t=%f <timegh>=%f sample=%d nh=%d\n",meas->Tnow,mean_gtime,meas->ngh,listh->usize);

	return;
}
/***************************************************
*  call routines that measure and store measures   *
*  during the time loop                            *
****************************************************/
void measures(TimeMeasures *meas){
	
        #ifdef DENSb1xT
	if((meas->Tnow>=meas->Ti)&&(meas->Tnow<=meas->Tf)){
		densB1Xt(meas);
	}
        #endif
        #ifdef AVERINVxT
	if((meas->Tnow>=meas->Ti)&&(meas->Tnow<=meas->Tf)){
		averInvestmentXt(meas);
	}
        #endif
        #ifdef MEANBFRACxT
	if((meas->Tnow>=meas->Ti)&&(meas->Tnow<=meas->Tf)){
		meanFracXt(meas);
	}
        #endif
        #ifdef SAVE_CONFIG
	double err=1e-7;
	if((meas->Tnow>=meas->saveT-err)&&(meas->Tnow<=meas->saveT+err)&&(meas->nfiles<NF)){
		printf("Time of measure:%f,  ",meas->Tnow);
		save_config(meas);
		meas->saveT=meas->Tnow+(double)NInterv*Dt_ref;
		printf("Next time:%f\n",meas->saveT);
		++meas->nfiles;
	}
        #endif
        #ifdef INV_DIST

	if((meas->Tnow>=500.)&&(meas->nfiles<NF)){
		printf("Time (measuring):%d\n",meas->NTnow);
		invDistXt(meas);
		++meas->nfiles;
	}else if(meas->Tnow<500.){
		printf("Time (before measuring starts):%d\n",meas->NTnow);
	}else{
		printf("Time (after measuring finishes):%d\n",meas->NTnow);
	}
        #endif
	#ifdef NUMHEVENTSxT
	numHostEventsPerDtXt(meas);
        #endif
	#ifdef CORRxT
	spatialCorrXt(meas);
        #endif
	#ifdef GENTIME
	genHostTime(meas);
        #endif


        return;
}

