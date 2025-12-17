/* meaures.c */

#include"globals.h"
#include"tools.h"
#include"measures.h"
/**********************************************
*   stores the frequency of beneficial        *
*   bacteria in a single host                 *
**********************************************/
void densB1Xt(double tnow,double ti,int idh){
	int nh;

	nh=listh->usize;

	if(tnow==ti){
		fprintf(fdensb1Xt,"#1:time 2:bac[single host,most helpful type] 3:dens_host \n");
	}
	
	fprintf(fdensb1Xt,"%f %f %f\n",tnow,bac[idh][TYPES-1]/spar->micr[idh],(double)nh/N);
       	fflush(fdensb1Xt);
       	printf("%f %f %f\n",tnow,bac[idh][TYPES-1]/spar->micr[idh],(double)nh/N);
        return;
}
/**********************************************
*   stores the frequency of beneficial       *
*   bacteria in a single host                *
**********************************************/
void averInvestmentXt(int ntnow,double tnow,int nti){
	int i,j,idh,nh;
	double averinv,averinvH,tot_micr;
	
	nh=listh->usize;

	if(ntnow==nti){
		fprintf(finvCumul,"#1:time 2:average cumulative investment 3:dens_host 4:number of time steps\n");
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
        fprintf(finvCumul,"%f %f %f %f %d\n",tnow,averinv,(double)tot_micr/nh,(double)nh/N,ntnow);
        fflush(finvCumul);
        printf("t=%f averinv=%f avmicrdens=%f nh/N=%f numsteps=%d\n",tnow,averinv,(double)tot_micr/nh,(double)nh/N,ntnow);

        return;
}
/**************************************
*  snapshot of host network           *
*  The colors indicate the acumulated *
*  investment of each host            *
***************************************/
void save_config(int tf,int tnow){
        int i,j,k,id;
        double pointsize,*w;
        char *name,*namedat,*name_gp;
        FILE *fconfig,*fgp;
	
        name=(char *)malloc(sizeof(char)*150);
        namedat=(char *)malloc(sizeof(char)*180);
        name_gp=(char *)malloc(sizeof(char)*180);

        w=(double *)calloc(N,sizeof(double));

        sprintf(name,"snapshot_N%d_Ty%d_Kh%d_Bv%f_Gh%d_cost%0.4f_sigma%0.2f_mu%f_sb%0.1f_mig%f_T%d",N,TYPES,spar->kh,Bacv,spar->gh,spar->cost,spar->sigma,spar->mu,spar->sb,spar->mig,tf);

        if(L>=100){
                pointsize=0.5;
        }else if(L>=50){
                pointsize=1.5;
        }else{
                pointsize=2;
        }
#if (SAVE_CONFIG_ID==0)
        sprintf(namedat,"%s_idt%f.dat",name,(double)tnow/tf);
        sprintf(name_gp,"%s_idt%f.gp",name,(double)tnow/tf);
#endif
#if (SAVE_CONFIG_ID==1)
        sprintf(namedat,"%s_Tf%d.dat",name,tnow);
        sprintf(name_gp,"%s_Tf%d.gp",name,tnow);
#endif

        fconfig = fopen(namedat,"w");
        fgp = fopen(name_gp,"w");


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

#if (SAVE_CONFIG_ID==0)
        fprintf(fgp,"set term png size 720,540\n");
       fprintf(fgp,"set output'%s_idt%f.png'\n",name,(double)tnow/tf);
#endif
#if (SAVE_CONFIG_ID==1)
        fprintf(fgp,"set term post eps enha color 20\n");
        fprintf(fgp,"set output'%s_Tf%d_fm.eps'\n",name,tnow);
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
	int i,j,nh;
	double avinvH;

	nh=listh->usize;
	for(i=0; i<nh; ++i){
		aveinvH=0.;
		for(j=0; j<TYPES; ++j){
			avinvH+=bac[i][j]*spar->inv[j]/spar->micr[j];
		}
		id=(int)(avinvH/binsize);
		++hist_inv[id]
	}


	return;
}
/**********************************************
*   this routine united all time measures     *
*   done inside the time loop                 *
**********************************************/
void measures(TimeMeasures meas){
	
	#ifdef TIME_VARS	
	int nh=listh->usize;
	for(int i=0; i<nh; ++i){
		int idh=listh->vec[i];
		double lowinv=0.;
		double medinv=0.;
		double highinv=0.;
		for(int j=0; j<TYPES;++j){
			if(spar->inv[j]<=0.1){
				lowinv+=bac[idh][j];
			}else if(spar->inv[j]<=0.5){
				medinv+=bac[idh][j];
			}else{
				highinv+=bac[idh][j];
			}
		}
		fprintf(fvarsXt,"%f %f %f %f %f\n",meas->Tnow,spar->micr[idh],lowinv,medinv,highinv);
		fflush(fvarsXt);
		printf("%f %f %f %f %f\n",meas->Tnow,spar->micr[idh],lowinv,medinv,highinv);
	}
	#endif
        #ifdef DENSb1xT
	if((meas.Tnow>=meas.Ti)&&(meas.NTnow<=meas.NTf)){
		densB1XtXt(meas.NTnow,meas.Ti,meas.idh_h);
	}
        #endif
        #ifdef AVERINVxT
	if((meas.NTnow>=meas.NTi)&&(meas.NTnow<=meas.NTf)){
		averInvestmentXt(meas.NTnow,meas.Tnow,meas.NTi);
	}
        #endif
        #ifdef SAVE_CONFIG
	if(meas.NTnow==meas.saveT){
		meas.saveT+=meas.interv;
		save_config(meas.NTnow,meas.NTf);

	}
        #endif


        return;
}

