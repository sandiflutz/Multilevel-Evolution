/*tools.c*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include"randgen_ufrgs.h"
#include"tools.h"

/********************************************************
 *                 Factorial                            *
 ********************************************************/
long int factorial(int n){
	int i,fac;
	
	if(n<=1){
		fac=1;
	}else{
		fac=n;
		for(i=n-1; i>0; --i){
			fac*=i;
		}
	}
		
	return fac;
}
/********************************************************
*   Normal Distribution: calculates and returns         *
*   the probability of x in a normal distribution       *
*   with standard deviation of 1 and mean of 0          *
********************************************************/
double normalProb(double x){
        double p;

        p=exp(-x*x/2.)/sqrt(2.*Pi);

        return p;
}
/****************************************************************
*      Abramowitz and Stegun approximation to 			*
*      the complementary cumulative distribution function	*
*      for a standart normal (and Horner's method). 		*
*      Source: 							*
*      https://www.johndcook.com/blog/normal_cdf_inverse/	*
*****************************************************************/
double rationalApprox(double t){
	double x,c[3],d[3];

	c[0] = 2.515517;
	c[1] = 0.802853;
	c[2] = 0.010328;

	d[0] = 1.432788;
	d[1] = 0.189269;
	d[2] = 0.001308;

	x = t-((c[2]*t+c[1])*t)+c[0]/(((d[2]*t + d[1])*t + d[0])*t + 1.);

	return x;
}
/****************************************************************
*     Return a value x from the inverse cumulative 		*
*     distribution function of a standart normal		*
*****************************************************************/
double invNormalCDF(double p){
	double t,finvcdf;

	if(p<0.){
		t=-sqrt(-2.*log(p));
	}else{
		t=sqrt(-2.*log(1-p));
	}
	
	finvcdf=rationalApprox(t);
	
	return finvcdf;
}
/****************************************************************
*     Return the cumulative distribution function for a normal 	*
*     distribution of x						*
*****************************************************************/
double normalCDF(double x){
	double fcdf;

	fcdf=(1.+erf(x/sqrt(2.)))/2.;

	return fcdf;
}
/****************************************************************
*     Generate a gaussian random number                         *
*****************************************************************/
double gaussRandNum(double mean, double sigma){
        double x,y,r,nr,fac;
               	
	do{
		x=2.*FRANDOM-1.;//random number between -1 and 1
		y=2.*FRANDOM-1.;
		r=x*x+y*y;
	}while(r>1. || r==0.);
                
	fac=sigma*sqrt(-2.*log(r)/r);
                
	x=fac*x+mean;
	y=fac*y+mean;

	if(FRANDOM<0.5){
 		nr=x;
	}else{
		nr=y;
	}
	
	return nr;
}
/******************************************************************
*     Generate a gaussian random number (from a truncated dist.)  *
*******************************************************************/
double truncGaussRandNum(double mean, double sigma,double xmin,double xmax){
        double nr;

        do{//reject samples outside of the range [xmin;xmax]
		nr=gaussRandNum(mean,sigma);
        } while( nr<xmin || nr > xmax);

        return nr;
}
/****************************************************************
*     	Generate a sample of normaly distributed		*
*     	random numbers in the range [xmin;xmax]			*
*****************************************************************/
void normalRandSample(double mean, double sigma,double xmin,double xmax,int sample,double *randvec){
	int numr;
        double x,y,r,fac;
               	
	numr=0;
        while(numr<sample){
		do{
			x=2.*FRANDOM-1.;//random number between -1 and 1
			y=2.*FRANDOM-1.;
			r=x*x+y*y;
		}while(r>1. || r==0.);
		fac=sigma*sqrt(-2.*log(r)/r);
                
		x=fac*x+mean;
		y=fac*y+mean;

		if((x>=xmin)&&(x<=xmax)){
			randvec[numr]=x;
			++numr;
		}
		if((numr<sample)&&(y>=xmin)&&(y<=xmax)){
			randvec[numr]=y;
			++numr;
		}
	}

        return;
}
/****************************************************************
*      Draw a random integer from a poisson distribution o      *
*      mean @lambda                                             *
*****************************************************************/
int poissonRandKnuth(double lambda){
	int k;
	double l,p,nr;

	k=0;
	p=1.;
	l=exp(-lambda);
	while(p>l){
		++k;
		nr=FRANDOM;
		p*=nr;
	}

	return k-1;
}
/************************************************************************************
*  Simple addition of an element to a list: add it at the end of the list and       *
*  increment the size of the list                                                   *
*************************************************************************************/
void listSimpleAdd(DynList *list,int add_elem){

	if((list->usize)<(list->size)){
		list->vec[list->usize]=add_elem;
		++list->usize;
	}else{
                exit(EXIT_FAILURE);
        }

	return;
}
/************************************************************************************
*  Simple subtraction of an element of a list: take the lest element of the         *
*  list and add it to the position of the element being subtracted and then         *
*  decrement the size of the list                                                   *
*************************************************************************************/
void listSimpleSub(DynList *list,int id_e){

	if(list->usize>0){
		list->vec[id_e]=list->vec[list->usize-1];
		--list->usize;
	}else{
                exit(EXIT_FAILURE);
        }

	return;
}
/***********************************************************************************
*      Add an element of a list (size=@size and currently used size=@usize).       *
*      Elements added to the part of the list being used come from the             *
*      second part of the list (elements from @usize to @size-1)                   * 
*      Def: @add_elem is the element being add to the list,@id_e is the old        *
*      position of this elemen#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>t (in the second part of the list)                   *  
************************************************************************************/
void listAdd(DynList *list,int add_elem,int id_e){
	int old_elem;

	if((list->usize)<(list->size)){
		old_elem=list->vec[list->usize];
		list->vec[list->usize]=add_elem;
		list->vec[id_e]=old_elem;
		++list->usize;
	}else{
                exit(EXIT_FAILURE);
        }

        return;
}
/***********************************************************************************
*      Subtract an element of a list (size=@size and currently used size=@usize).  *
*      Elements subtracted from the list are stored in the second part of the list *
*      (elements from @usize to @size-1)                                           *
*      Def: @sub_elem is the element being subtracted from the list, @id_e is the  *
*      old position of this element (in the first part of the list)                *
************************************************************************************/
void listSub(DynList *list,int sub_elem,int id_e){
	int old_elem;

	if(list->usize>0){
		old_elem=list->vec[id_e];
		list->vec[id_e]=list->vec[list->usize-1];
		list->vec[list->usize-1]=old_elem;
		--list->usize;
	}else{
                exit(EXIT_FAILURE);
        }

        return;
}
/****************************************************************
*	exchange 2 elements of a list of integers		*
*****************************************************************/
void exchange(int *list, int id1, int id2){
	int old;

	old=list[id1];
	list[id1]=list[id2];
	list[id2]=old;

	return;
}
/****************************************************************
*	exchange 2 elements of a list of doubles		*
*****************************************************************/
void exchangeF(double *list, int id1, int id2){
	double old;

	old=list[id1];
	list[id1]=list[id2];
	list[id2]=old;

	return;
}
/************************************************************************************
*    find the number of live neighbors of a specific host using the neighbors       *
*    matrix and the list of live hosts. Store live neighbors positions              *
************************************************************************************/
void searchLiveNeighbors(int sa,int idh,int *state,int **neighbors,DynList *alive_viz){
	int j,idv,ne,viz,nviz;

	nviz=0;
	ne=0;

	viz=alive_viz->size-1;//the full size of the list includes the focus site 
	for(j=0;j<viz;++j){
		idv=neighbors[idh][j];
		if(state[idv]==sa){//sa is the state alive state
			alive_viz->vec[nviz]=idv;
			++nviz;
		}else{
			alive_viz->vec[viz-1-ne]=idv;//storing empty neighboring sites at the end of each "line"
			++ne;
		}
	}
	alive_viz->vec[nviz]=idh;
	alive_viz->usize=nviz+1;


	return;
}
/************************************************************************************
*    find the number of empty sites in the neighborhood of a specific site using    * 
*    the neighbors matrix. Store the neighboring empty sites positions              *
************************************************************************************/
void searchEmptyNeighbors(int se,int idh,int *state,int **neighbors,DynList *empty_viz){
	int j,idv,ne,viz;

	ne=0;
	viz=empty_viz->size;
	for(j=0; j<viz; ++j){
		idv=neighbors[idh][j];
		if(state[idv]==se){//se is the dead state
			empty_viz->vec[ne]=idv;
			++ne;
		}
	}
	empty_viz->usize=ne;
	
	return;
}
/************************************************************************************
*    find the number of live neighbors using the neighbors matrix and              *
*    the list of live hosts. Store live neighbors positions                        *
*    If the nh*maxviz<nh/2 choose this instead (where nh is the number if live     * 
*    hosts ans maxviz in the number of neighbors)                                   *
************************************************************************************/
void searchAllLiveNeighbors(int dead,int sites,int nh,int maxviz,int *listh,int *host,int **neighbor,int *nviz,int *list_aviz){
	int i,j,idh,idv,idv_e,ne;

	memset(nviz,0,sizeof(int)*sites);

	for(i=0; i<nh; ++i){
		idh=listh[i];
		ne=0;
		for(j=0;j<maxviz-1;++j){//position maxviz-1 is occupied by the if of the focus host
			idv=neighbor[idh][j];
			if(host[idv]!=dead){
				list_aviz[idh*sites+nviz[idh]]=idv;
				++nviz[idh];
			}else{
				idv_e=maxviz-1-ne;
				list_aviz[idh*sites+idv_e]=idv;//storing empty neighboring sites at the end of each "line"
				++ne;
			}
		}
		list_aviz[idh*sites+nviz[idh]]=idh;
		++nviz[idh];
	}

	return;
}
/************************************************************************************
*    find the number of live neighbors using the links matrix and                   *
*    the list of live hosts. Store live neighbors positions.                        *
*    If the nh*maxviz>nh/2 choose this instead (where nh is the number if live      * 
*    hosts ans maxviz in the number of neighbors)                                   *
*************************************************************************************/
void searchAllLiveLinks(int sites,int nh,int *listh,int *netlink,int *nviz,int *list_aviz){
	int i,j,idh1,idh2;

	memset(nviz,0,sizeof(int)*sites);
	for(i=0; i<nh-1; ++i){
		idh1=listh[i];
		for(j=i+1; j<nh; ++j){
			idh2=listh[j];
			if(netlink[idh1*sites+idh2]==1){
				list_aviz[idh1*sites+nviz[idh1]]=idh2;
				++nviz[idh1];
				list_aviz[idh2*sites+nviz[idh2]]=idh1;
				++nviz[idh2];
			}
		}
	}

	return;
}
/***************************************************
* Set neighbors for a square lattice network.      *
* Indexes order: neighbors[site_id][k-th neighbor] *
****************************************************/
void squareLattice(int **neighbor,int viz,int lsize){
        int i,k,sites,base;

        /*neighbors order:
         * 0=top,1=right,2=botton,3=left;
         * Diagonals:4=left-top,5=right-top,6=right-bottom,7=left-bottom */

	sites=lsize*lsize;
	k=0;
        while(k<viz){
                for(i=0; i<sites; ++i){
                        base=(i/lsize)*lsize;
                        switch(k){
                                case 0: neighbor[i][k] = (i-lsize+sites)%sites;//up
                                        break;
                                case 1: neighbor[i][k]= base + (i+1)%lsize;//right
                                        break;
                                case 2: neighbor[i][k]= (i+lsize)%sites;//down
                                        break;
                                case 3: neighbor[i][k]= base + (i-1+lsize)%lsize;//left
                                        break;
                                case 4: neighbor[i][k]= base + (i-1+lsize)%lsize;//(left) left-up
                                        neighbor[i][k]= (neighbor[i][k]-lsize+sites)%sites;//(up)
                                        break;
                                case 5: neighbor[i][k]= base + (i+1)%lsize;//(right) right-up
                                        neighbor[i][k]=(neighbor[i][k]-lsize+sites)%sites;//(up)
                                        break;
                                case 6: neighbor[i][k]= base + (i+1)%lsize;//(right) right-down
                                        neighbor[i][k]= (neighbor[i][k]+lsize)%sites;//(down)
                                        break;
                                case 7: neighbor[i][k]= base + (i-1+lsize)%lsize;//(left)  left-down
                                        neighbor[i][k]= (neighbor[i][k]+lsize)%sites;//(down)
                                        break;
                        }
                }
                ++k;
        }

        return;
}
/**********************************************************
*               Set network of a complete graph           *
*               (without self-neighboring)
***********************************************************/
void setCompleteGraph(int **neighbor,int sites){
        int i,j,*nviz;

        nviz=(int *)calloc(sites,sizeof(int));
        memset(nviz,0,sizeof(int)*sites);

        for(i=0; i<sites-1; ++i){
                for(j=i+1; j<sites; ++j){
                        neighbor[i][nviz[i]]=j;
                        ++nviz[i];
                        neighbor[j][nviz[j]]=i;
                        ++nviz[j];
                }
        }

        free(nviz);
        return;
}
/***********************************************************
*       Set a neighbors classification vector:             *
*       clneighbor[i][j]=k: j is the k-th neighbor of i    *
************************************************************/
void classifyNeighbors(int **neighbor,int **clneighbor,int nviz,int sites){
	int i,k;

	if(neighbor!=NULL){
		for(i=0; i<sites; ++i){
			for(k=0; k<nviz; ++k){
				clneighbor[i][neighbor[i][k]]=k;
			}
		}
	}else{
		 exit(1);
	}
	

        return;
}
/**************************************************
*       Randomly select a neighbor from           *
*       a neighbors list, from indexes 1d1 to id2 *
***************************************************/
int randNeighborID(int vec_id,int *vec,int id1,int id2,int size){
	int id_viz,idmax,idmin;

	if(id1<id2){
		idmin=id1;
		if(id2>size)id2=size;
		idmax=id2;
	}else if(id2<id1){
		idmin=id2;
		if(id1>size)id1=size;
		idmax=id1;
	}else{
		idmin=id1;
		idmax=idmin;
	}

	do{
		id_viz=(int)(FRANDOM*(idmax-idmin))+idmin;
	}while(vec[id_viz]==vec_id);

	return id_viz;
}
/*************************************************
*       Set network links between sites          *
**************************************************/
void setNetLinks(int **neighbor,int *netlink,int nviz,int sites){
        int i,j,viz;

        memset(netlink,0,sizeof(int)*sites*sites);
        for(i=0; i<sites; ++i){
                for(j=0; j<nviz; ++j){
                        viz=neighbor[i][j];
                        netlink[i*sites+viz]=1;
                }
        }

        return;
}
/****************************************************
* Calculates the cumulative probability vector from *
* a simple probability vector                       *
*****************************************************/
void cumulProb(int size, double *prob,double *cprob){
        int i;

        cprob[0]=prob[0];
        for(i=1; i<size; ++i){
                cprob[i]=cprob[i-1]+prob[i];
        }

        return;
}
/*********************************************************
 *   return the smallest index k that satisfy vec[k]>nr  *                           
 ********************************************************/
int upperBoundStrict(double nr,double *vec, int vsize){
	int left, right,mid;
		
	left=0;
	right=vsize;
	while(left<right){
		mid=left+(right-left)/2;
		if(vec[mid]>nr){
			right=mid;
		}else{
			left=mid+1;
		}
	}

	if(left>=vsize)left=vsize-1;

	return left;
}
/*********************************************************
*         randomly selects and event (towering method):  *
*   1)calculates the cumulative probabilitities          *
*   2)Randomly selects a number between 0 and the last   *
*   element of the cumulative probability vector         *
*   3)calls the bissectionSearch routine to find the     *
*   position, in the vector, the number belongs to.      *
*********************************************************/
int selectEvent(double nr,double *event,int size){
        int id;
        double *cprob;

	cprob=(double *)calloc(size,sizeof(double));
	cumulProb(size,event,cprob);

      //  nr=FRANDOM*cprob[size-1];
        id=upperBoundStrict(nr, cprob, size);

	free(cprob);
        return id;
}
/***********************************************
*         randomly selects and event           *
*         using the respective cumulative      *
*         probability vector                   *
************************************************/
int selectEventCP(double nr,double *cprob,int size){
        int id;

      //  nr=FRANDOM*cprob[size-1];
        id=upperBoundStrict(nr, cprob, size);

	return id;
}
/***********************************************
*           Find the largest element           *
*           of a vector and returns            *
*           its position in the vector         *
************************************************/
int findMaxElement(int size,double *max_e,double *vec){
        int i,idm;

        *max_e=vec[0];
	idm=0;
        for(i=1; i<size; ++i){
                if(vec[i]>*max_e){
			*max_e=vec[i];
			idm=i;
		}
        }

        return idm;
}
/********************************************************
 * return the smallest number between the @a and @b     *
 ********************************************************/
double minNum(double a, double b){
	double min;

	if(a<=b){
		min=a;
	}else{
		min=b;
	}


	return min;
}
/****************************************************
 * Create fill a vector with log spaced elements    *
 * 	- range: [m0*10^(e0),mf*10^(ef)]            *
 * 	- number of elements: size                  * 
 ****************************************************/
void logSpacedVec(double *vec,double e0,double ef,double m0, double mf, int size){
	int i;
	double step=(ef+log10(mf)-e0-log10(m0))/(size-1);

	for(i=0; i<size; ++i){
		vec[i]=pow(10.,(e0+log10(m0)+i*step));
	}

	return;
}
/********************************************************
 *              Build a histogram                       *
 ********************************************************/
void buildHistogram(double binsize,int nbins,int vsize,double *vec,int *hist){
        int i,id;

        for(i=0; i<nbins; ++i){
                hist[i]=0;
        }
        for(i=0; i<vsize; ++i){
                id=(int)(vec[i]/binsize);
                ++hist[id];
        }

        return;
}
/********************************************************
*  calculate the root mean square of the elements	* 
*  from index idi to index idf of a vector vec[]	*
*********************************************************/
void calcRMSError(double *vec,int idi, int idf,double *stats){
	int i,sample;
	double aver=0., err=0.;

	if(!stats)exit(1);

	sample=idf-idi;
	for(i=idi; i<=idf; ++i){
		aver+=vec[i];
	}
	aver/=(double)sample;
	for(i=idi; i<=idf; ++i){
		err+=(vec[i]-aver)*(vec[i]-aver);
	}
	err=sqrt(err/(double)sample);
	
	stats[0]=aver;
	stats[1]=err;

	return;
}
/***************************************************
*  calculate spatial corretation for a square      *
*  lattice type of system                          *
****************************************************/
double spatialCorr(int *state,int sites, int dist,int right,int down,int **neighbor,double *corr){
        int i,j,idvx,idvy;
        double corrx=0.;
        double corry=0.;
        double corr_tot=0.;

        for(i=0; i<sites; ++i){
                idvx=i;
                idvy=i;
                for(j=0; j<dist; ++j){
                        idvx=neighbor[idvx][right];
                        idvy=neighbor[idvy][down];
                }
                corrx+=(double)state[i]*state[idvx];
                corry+=(double)state[i]*state[idvy];
        }
	corr_tot=(corrx+corry)/(2.*sites);
	corrx/=(double)sites;
	corry/=(double)sites;
	
	corr[0]=corrx;
	corr[1]=corry;

        return corr_tot;
}
/***************************************************
*  calculate spatial corretation for a square      *
*  lattice type of system in 1 direction only      *
****************************************************/
double spatialCorr1d(int *state,int sites, int dist,int id_direction,int **neighbor){
        int i,j,idv;
        double corr_dir=0.;

        for(i=0; i<sites; ++i){
                idv=i;
                for(j=0; j<dist; ++j){
                        idv=neighbor[idv][id_direction];
                }
                corr_dir+=(double)state[i]*state[idv];
        }
	corr_dir/=(double)sites;

        return corr_dir;
}
/****************************************
 * 	Heaviside Step Function		*
 * 	for x in (-1.,1]		*
 ***************************************/
double stepFuncBounded(double x){
	return (ceil(x));
}
