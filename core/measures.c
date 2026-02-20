/* meaures.c */

#include"globals.h"
#include"tools.h"
#include"measures.h"
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

	if(meas->NTnow==meas->NTi){
		fprintf(meas->file_tmeas,"#1:time 2:average investment level on an isolated host\n");
	}
	
	averinv=0.;
	for(j=0; j<TYPES; ++j){
		averinv+=bac[meas->idh_h1][j]*spar->s[j]*spar->inv[j]/spar->micr[meas->idh_h1];
	}
	
	fprintf(meas->file_tmeas,"%f %f %f %f %f\n",meas->Tnow,averinv,bac[meas->idh_h1][0]/spar->micr[meas->idh_h1],bac[meas->idh_h1][1]/spar->micr[meas->idh_h1],spar->micr[meas->idh_h1]);
	printf("time=%f averinv=%f fbac0=%f fbac1=%f micr=%f\n",meas->Tnow,averinv,bac[meas->idh_h1][0]/spar->micr[meas->idh_h1],bac[meas->idh_h1][1]/spar->micr[meas->idh_h1],spar->micr[meas->idh_h1]);
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
	
	if(meas->NTnow==meas->NTi){
		fprintf(meas->file_tmeas,"#1:time 2:mean freq. of bac of type=%d ",TYPES-1);
		for(j=TYPES-2; j>=0; --j){
			fprintf(meas->file_tmeas,"%d:mean freq. of bac of type=%d ",TYPES-j+1,j);
		}
		fprintf(meas->file_tmeas,"\n");
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
        fprintf(meas->file_tmeas,"\n");
        printf("\n");

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
	
	if(meas->NTnow==meas->NTi){
		fprintf(meas->file_tmeas,"#1:time 2:average cumulative investment 3:dens_host 4:number of time steps\n");
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
        int i,j,id,idlisth,nh,namelen,dnl;
        double pointsize,*freqInvH;
        char *name,*namedat,*name_gp,*nt_format;
        FILE *fconfig,*fgp;
	
	/****creating files (if there are file with the same names, there are subscribed)********/
	dnl=100;
        namelen=strlen(meas->ftname_pars)+strlen(fdatapath)+dnl;
        nt_format=(char *)malloc(sizeof(char)*dnl);//format of the time term, on the name (files for animations have to be ordered: better to use a 0<t<1 format)
        name=(char *)malloc(sizeof(char)*namelen);//base name
        namedat=(char *)malloc(sizeof(char)*namelen+dnl);//data file name
        name_gp=(char *)malloc(sizeof(char)*namelen+dnl);//name for the gnuplot data file

        sprintf(name,"%ssnapshot_%s",fdatapath,meas->ftname_pars);

        if(L>=100){
                pointsize=0.5;
        }else if(L>=50){
                pointsize=1.5;
        }else{
                pointsize=2;
        }
#if (FIG_EXT==0)
        sprintf(nt_format,"idt%f",(double)meas->NTnow/meas->NTf);
#elif
        sprintf(nt_format,"Tf%d",meas->NTnow);
#endif
        sprintf(namedat,"%s%s_%s.dat",fdatapath,name,nt_format);
        sprintf(name_gp,"%s%s_%s.gp",fdatapath,name,nt_format);

        fconfig = fopen(namedat,"w");
        fgp = fopen(name_gp,"w");
	

	/*******setting average investment per host vector and filling data file***********************/

	nh=listh->usize;
	freqInvH=(double *)calloc(nh,sizeof(double));/*indexes are the same as in @listh->vec 
						      * (real host indexes are idh=listh->vec[i] (for 0<=i<nh)) */
	calcInvFreq(freqInvH);
	
	idlisth=0;
        for(i = 0; i < L; ++i){//linha
                for(j = 0; j < L; ++j){//coluna
                        id=j+i*L;
			if(id<listh->vec[idlisth]){//@id not on the list of alive hosts
                                fprintf(fconfig,"%d %d %d\n",i,j,-1);
			}else{//id==listh->vec[idlisth], @id is alive
                                fprintf(fconfig,"%d %d %f\n",i,j,freqInvH[idlisth]);
				++idlisth;
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
	fprintf(fgp,"set output'%s_idt%f.png'\n",name,(double)meas->NTnow/meas->NTf);
#elif
        fprintf(fgp,"set term post eps enha color 20\n");
        fprintf(fgp,"set output'%s_Tf%d_fm.eps'\n",name,meas->NTnow);
#endif

        fprintf(fgp,"unset key\n");
        fprintf(fgp,"unset xtics\n");
        fprintf(fgp,"unset ytics\n");
        fprintf(fgp,"L=%d\n",L);
        fprintf(fgp,"set yr[0:L-1]\n");
        fprintf(fgp,"set xr[0:L-1]\n");
        fprintf(fgp,"set size square\n");
        fprintf(fgp,"set pointsize %0.2lf\n",pointsize);
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
        nt_format=(char *)malloc(sizeof(char)*dnl);//base name
        name=(char *)malloc(sizeof(char)*namelen);//format of the time term, on the name (files for animations have to be ordered: better to use a 0<t<1 format)
        namedat=(char *)malloc(sizeof(char)*namelen+dnl);//data file name
        name_gp=(char *)malloc(sizeof(char)*namelen+dnl);//name for the gnuplot data file
        
	sprintf(name,"invDistXt_%s",meas->ftname_pars);

#if (FIG_EXT==0)
        sprintf(nt_format,"idt%f",(double)meas->NTnow/meas->NTf);
#elif
        sprintf(nt_format,"Tf%d",meas->NTnow);
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
	
	if(meas->NTnow==meas->NTi)fprintf(fhist,"#1:average inv. in hosts 2:frac. of hosts 3:numsteps\n");
	for(i=0; i<nbins; ++i){
		fprintf(fhist,"%f %f %f %d\n",(double)i*binsize,(double)hist_inv[i]/nh,freqInvH[0],meas->NTnow);
	}
	
	/*******************filling gnuplot file*****************************/

#if (FIG_EXT==0)
       fprintf(fgp,"set term png size 720,540\n");
       fprintf(fgp,"set output'%s_idt%f.png'\n",name,(double)meas->NTnow/meas->NTf);
#elif
        fprintf(fgp,"set term post eps enha color 20\n");
        fprintf(fgp,"set output'%s_Tf%d_fm.eps'\n",name,meas->NTnow);
#endif

	fprintf(fgp,"set title'{/=15 time steps=%d}'\n",meas->NTnow);
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
*  call routines that measure and store measures   *
*  during the time loop                            *
****************************************************/
void measures(TimeMeasures *meas){
	
        #ifdef DENSb1xT
	if((meas->NTnow>=meas->NTi)&&(meas->NTnow<=meas->NTf)){
		densB1Xt(meas);
	}
        #endif
        #ifdef AVERINVxT
	if((meas->NTnow>=meas->NTi)&&(meas->NTnow<=meas->NTf)){
		averInvestmentXt(meas);
	}
        #endif
        #ifdef MEANBFRACxT
	if((meas->NTnow>=meas->NTi)&&(meas->NTnow<=meas->NTf)){
		meanFracXt(meas);
	}
        #endif
        #ifdef SAVE_CONFIG
	if((meas->NTnow>=meas->NTi)&&(meas->nfiles<NF)){
		printf("Time (measuring):%d\n",meas->NTnow);
		save_config(meas);
		++meas->nfiles;
	}
	printf("Time (before measuring starts):%d\n",meas->NTnow);
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


        return;
}

