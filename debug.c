/*debug.c*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"globals.h"
#include"randgen_ufrgs.h"
#include"tools.h"
#include"evo.h"

int main(void){

	(void)start_randomic(0);

#ifdef DEBUG_CALC_DNT
	int dnt_h,dnt_b,dnt;
	double dth,dtb;

	dtb=0.05;
	dth=dtb/2.5;
	calcNumSteps(&dnt_h,&dnt_b,&dnt,dth,dtb);
	printf("dtb=%f dth=%f dnt_h=%d dnt_b=%d dnt=%d\n",dtb,dth,dnt_h,dnt_b,dnt);
	dth=dtb/1.1;
	calcNumSteps(&dnt_h,&dnt_b,&dnt,dth,dtb);
	printf("dtb=%f dth=%f dnt_h=%d dnt_b=%d dnt=%d\n",dtb,dth,dnt_h,dnt_b,dnt);
	dth=dtb/1.6;
	calcNumSteps(&dnt_h,&dnt_b,&dnt,dth,dtb);
	printf("dtb=%f dth=%f dnt_h=%d dnt_b=%d dnt=%d\n",dtb,dth,dnt_h,dnt_b,dnt);
#endif
#ifdef DEBUG_TRUNCNORM
	int i;
	double mean=0.01;
	double var=0.05;
	int bsize=0.01;
	int sample=10000;
	int beans=100;
	bsize=1./(double)beans;
	int *hist=(int *)calloc(beans,sizeof(int));
	memset(hist,0,sizeof(int)*beans);

	printf("#mean=%f  var=%f\n",mean,var);
	for(i=0; i<sample; ++i){
		double nr=gaussRandNum(mean,var,0.,1.);
		int id=(int)(nr/bsize);
		++hist[id];
		
	}

	for(i=0; i<beans; ++i){
		printf("%f %d %f\n",(double)i*bsize,i,(double)hist[i]/sample);
	}

	free(hist);
#endif
#ifdef DEBUG_POISSON
	int i,k,sample=1000;
	double lambda=1e-03;

	
	for(i=0; i<sample; ++i){
		k=poissonRandNum(lambda);
		printf("i=%d k=%d\n",i,k);
	}
	
#endif
#ifdef DEBUG_BISS
	int i,id;
	int size=10;
	int sample=1000;
	double nr;
	
	int *num=(int *)calloc(size,sizeof(int));
	double *vec=(double *)calloc(size,sizeof(double));
	double *cprob=(double *)calloc(size,sizeof(double));
	memset(vec,0.,sizeof(double)*size);
	memset(num,0.,sizeof(int)*size);
	vec[2]=0.2;
	vec[4]=0.2;
	vec[7]=0.1;
	vec[9]=0.5;
	cprob[0]=vec[0];
	for(i=1; i<10; ++i){
		cprob[i]=cprob[i-1]+vec[i];	
	}
	
	for(i=0; i<sample; ++i){
		nr=FRANDOM;
		id=bissectionSearch(nr,cprob,size);
		++num[id];
	}
	for(i=0; i<size; ++i){
		printf("i:%d prob[i]:%f vec[i]:%f cprob[i]:%f\n",i,(double)num[i]/sample,vec[i],cprob[i]);
	}
	free(vec);
	free(num);
	free(cprob);
#endif
	return 0;
}
