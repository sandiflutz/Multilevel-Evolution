/* meaures.c */

#include"globals.h"
#include"tools.h"
#include"measures.h"
/**********************************************
*   stores the frequency of beneficial        *
*   bacteria in a single host                 *
**********************************************/
void densB1Xt(TimeMeasures *meas){
	int nh;

	nh=listh->usize;

	if(meas->NTnow==meas->NTi){
		fprintf(meas->file_tmeas,"#1:time 2:bac[single host,most helpful type] 3:dens_host \n");
	}
	
	fprintf(meas->file_tmeas,"%f %f %f\n",meas->Tnow,bac[meas->idh_h1][TYPES-1]/spar->micr[meas->idh_h1],(double)nh/N);
       	printf("%f %f %f\n",meas->Tnow,bac[meas->idh_h1][TYPES-1]/spar->micr[meas->idh_h1],(double)nh/N);
        return;
}
/**********************************************
*   stores the frequency of beneficial       *
*   bacteria in a single host                *
**********************************************/
void averInvestmentXt(TimeMeasures *meas){
	int i,j,idh,nh;
	double averinv,averinvH,tot_micr;
	
	nh=listh->usize;

	if(meas->NTnow==meas->NTi){
		fprintf(meas->file_tmeas,"#1:time 2:average cumulative investment 3:dens_host 4:number of time steps\n");
	}

	averinv=0.;
	tot_micr=0.;
        for(i=0; i<nh; ++i){
                idh=listh->vec[i];

		tot_micr+=spar->micr[idh];
                averinvH=0.;
                for(j=0; j<TYPES; ++j){
                        averinvH+=spar->inv[j]*bac[idh][j];//the investment for strategy j=0 (neutral bacteria) is 0, so only helpers are being acounted for here
                }
                averinv+=averinvH;
        }
        averinv/=tot_micr;
        fprintf(meas->file_tmeas,"%f %f %f %f %d\n",meas->Tnow,averinv,(double)tot_micr/nh,(double)nh/N,meas->NTnow);
        printf("t=%f averinv=%f avmicrdens=%f nh/N=%f numsteps=%d\n",meas->Tnow,averinv,(double)tot_micr/nh,(double)nh/N,meas->NTnow);

        return;
}
/**************************************
*  snapshot of host network           *
*  The colors indicate the acumulated *
*  investment of each host            *
***************************************/
void save_config(TimeMeasures *meas){
        int i,j,k,id,nh,namelen,dnl;
        double pointsize,*w;
        char *name,*namedat,*name_gp;
        FILE *fconfig,*fgp;
	
	/****creating files (if there are file with the same names, there are subscribed)********/
	dnl=100;
        namelen=strlen(meas->ftname_pars)+strlen(fdatapath)+dnl;
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
        sprintf(namedat,"%s_idt%f.dat",name,(double)meas->NTnow/meas->NTf);
        sprintf(name_gp,"%s_idt%f.gp",name,(double)meas->NTnow/meas->NTf);
#elif
        sprintf(namedat,"%s_Tf%d.dat",name,meas->NTnow);
        sprintf(name_gp,"%s_Tf%d.gp",name,meas->NTnow);
#endif

        fconfig = fopen(namedat,"w");
        fgp = fopen(name_gp,"w");
	

	/*******setting average investment per host vector and filling data file***********************/

	nh=listh->usize;
        w=(double *)calloc(nh,sizeof(double));

        for(i = 0; i < L; ++i){//linha
                for(j = 0; j < L; ++j){//coluna
                        id=j+i*L;
                        w[id]=-1.;
                        if(spar->micr[id]>0.){
                                w[id]=0.;
                                for(k=0; k<TYPES; ++k){
                                        w[id]+=spar->inv[k]*spar->s[k]*bac[id][k]/spar->micr[id];
                                }
                        }

                        if(spar->micr[id]>0.){
                                fprintf(fconfig,"%d %d %f %f\n",i,j,w[id],spar->inv[TYPES-1]*bac[id][TYPES-1]/spar->micr[id]);
                        }else{
                                fprintf(fconfig,"%d %d %f %f\n",i,j,w[id],-1.);
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
        fprintf(fgp,"plot '%s' u 1:(($3>=0)?$2:1/0):(rgb(255*(1-$3),0,255*$3)) w p pt 5 lc rgb variable,\\\n",namedat);
        fprintf(fgp,"'%s' u 1:(($3<0)?$2:1/0):(rgb(0,0,0)) w p pt 5 lc rgb variable\n",namedat);


	/***freeing allocated memory and closing files******/
        fclose(fconfig);
        fclose(fgp);
        free(name_gp);
        free(namedat);
        free(name);
        free(w);

        return;
}
/**********************************************
 * calculates investment distribution among   *
 * hosts                                      *
 **********************************************/
void calcInvDist(double *hist_inv,double binsize){
	int i,idh,j,nh,id;
	double avinvH;

	nh=listh->usize;
	for(i=0; i<nh; ++i){
		idh=listh->vec[i];
		avinvH=0.;
		for(j=0; j<TYPES; ++j){
			avinvH+=bac[idh][j]*spar->inv[j]/spar->micr[idh];
		}
		id=(int)(avinvH/binsize);
		++hist_inv[id];
	}


	return;
}
/*****************************************************
 * store investment distribution among hosts and     *
 * create a gnuplot script to create graphics        *
 *****************************************************/
void invDistXt(TimeMeasures *meas){
	int i,nh,nbins=100,namelen,dnl;
	double binsize,*hist_inv;
        char *name,*namedat,*name_gp;
        FILE *fhist,*fgp;
	
	/****creating files (if there are file with the same names, there are subscribed)********/
	dnl=100;
	namelen=strlen(meas->ftname_pars)+strlen(fdatapath)+dnl;
        name=(char *)malloc(sizeof(char)*namelen);//base name
        namedat=(char *)malloc(sizeof(char)*namelen+dnl);//data file name
        name_gp=(char *)malloc(sizeof(char)*namelen+dnl);//name for the gnuplot data file
        
	sprintf(name,"%sinvDistXt_%s",fdatapath,meas->ftname_pars);

#if (FIG_EXT==0)
        sprintf(namedat,"%s_idt%f.dat",name,(double)meas->NTnow/meas->NTf);
        sprintf(name_gp,"%s_idt%f.gp",name,(double)meas->NTnow/meas->NTf);
#elif
        sprintf(namedat,"%s_Tf%d.dat",name,meas->NTnow);
        sprintf(name_gp,"%s_Tf%d.gp",name,meas->NTnow);
#endif
        fhist = fopen(namedat,"w");
        fgp = fopen(name_gp,"w");
	
	/*******setting histogram vector***********************/
	
	hist_inv=(double *)calloc(nbins,sizeof(double));
	memset(hist_inv,0.,sizeof(double)*nbins);
	binsize=1./(double)nbins;

	calcInvDist(hist_inv,binsize);

	/*******************filling data file*****************************/
	nh=listh->usize;

	if(meas->NTnow==meas->NTi)fprintf(meas->file_tmeas,"#1:time 2:hist_inv[type] 3:type 4:numsteps\n");
	for(i=0; i<nbins; ++i){
		fprintf(fhist,"%f %f %d\n",(double)i*binsize,(double)hist_inv[i]/nh,meas->NTnow);
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
	fprintf(fgp,"plot '%s' u 1:2 w l lw 3 lc 1\n",namedat);

	/**freeing allocated memory and closing files**/

	free(hist_inv);
	free(name);
	free(namedat);
	free(name_gp);
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
	if((meas->NTnow>=meas->Ti)&&(meas->NTnow<=meas->NTf)){
		densB1XtXt(meas);
	}
        #endif
        #ifdef AVERINVxT
	if((meas->NTnow>=meas->NTi)&&(meas->NTnow<=meas->NTf)){
		averInvestmentXt(meas);
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

	if((meas->NTnow>=meas->NTi)&&(meas->nfiles<NF)){
		printf("Time (measuring):%d\n",meas->NTnow);
		invDistXt(meas);
		++meas->nfiles;
	}else{
		printf("Time (before measuring starts):%d\n",meas->NTnow);
	}
        #endif


        return;
}

