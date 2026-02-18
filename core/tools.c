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
*     Generate a gaussian random number                         *
*****************************************************************/
double gaussRandNum(double mean, double var){
        double x,y,r,nr,fac;

               	
	do{
		x=2.*FRANDOM-1.;//random number between -1 and 1
		y=2.*FRANDOM-1.;
		r=x*x+y*y;
                
	}while(r>1. || r==0.);
                
	fac=var*sqrt(-2.*log(r)/r);
                
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
double truncGaussRandNum(double mean, double var,double a,double b){
        double x,y,r,nr,fac;

        do{//reject samples outside of the range [a;b]
                do{
                        x=2.*FRANDOM-1.;//random number between -1 and 1
                        y=2.*FRANDOM-1.;

                        r=x*x+y*y;
                }while(r>1. || r==0.);

                fac=var*sqrt(-2.*log(r)/r);

                x=fac*x+mean;
                y=fac*y+mean;

                if(FRANDOM<0.5){
                        nr=x;
                }else{
                        nr=y;
                }

        } while( nr<=a || nr >= b);

        return nr;
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
/**************Double List Version***********************************/
void listSimpleAddF(DynListF *list,double add_elem){

	if((list->usizef)<(list->sizef)){
		list->vecf[list->usizef]=add_elem;
		++list->usizef;
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
/**************Double List Version***********************************/
void listSimpleSubF(DynListF *list,int id_e){

	if(list->usizef>0){
		list->vecf[id_e]=list->vecf[list->usizef-1];
		--list->usizef;
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
/**************Double List Version***********************************/
void listAddF(DynListF *list,double add_elem,int id_e){
	double old_elem;

	if((list->usizef)<(list->sizef)){
		old_elem=list->vecf[list->usizef];
		list->vecf[list->usizef]=add_elem;
		list->vecf[id_e]=old_elem;
		++list->usizef;
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
/**************Double List Version***********************************/
void listSubF(DynListF *list,double sub_elem,int id_e){
	double old_elem;

	if(list->usizef>0){
		old_elem=list->vecf[id_e];
		list->vecf[id_e]=list->vecf[list->usizef-1];
		list->vecf[list->usizef-1]=old_elem;
		--list->usizef;
	}else{
                exit(EXIT_FAILURE);
        }
	

        return;
}
/************************************************************
*  		exchange 2 elements of a list               *
*************************************************************/
void exchange(int *list, int id1, int id2){
	int old;

	old=list[id1];
	list[id1]=list[id2];
	list[id2]=old;

	return;
}
/************************************************************************************
*    find the number of live neighbors of a specific host using the neighbors       *
*    matrix and the list of live hosts. Store live neighbors positions              *
************************************************************************************/
void searchLiveNeighbors(int dead,int idh,int *host,int **neighbor,DynList *alive_viz){
	int j,idv,ne,viz,nviz;

	nviz=0;
	ne=0;

	viz=alive_viz->size-1;//the full size of the list includes the focus site 
	for(j=0;j<viz;++j){
		idv=neighbor[idh][j];
		if(host[idv]!=dead){
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
void squareLattice(int **neighbor,int viz,int sites){
        int i,k=0,l,base;

        l=sites/2;

        /*neighbors order:
         * 0=top,1=right,2=botton,3=left;
         * Diagonals:4=left-top,5=right-top,6=right-bottom,7=left-bottom */

        while(k<viz){
                for(i=0; i<sites; ++i){
                        base=(i/l)*l;
                        switch(k){
                                case 0: neighbor[i][k]= (i-l+sites)%sites;//up
                                        break;
                                case 1: neighbor[i][k]= base + (i+1)%l;//right
                                        break;
                                case 2: neighbor[i][k]= (i+l)%sites;//down
                                        break;
                                case 3: neighbor[i][k]= base + (i-1+l)%l;//left
                                        break;
                                case 4: neighbor[i][k]= base + (i-1+l)%l;//(left) left-up
                                        neighbor[i][k]= (neighbor[k][i]-l+sites)%sites;//(up)
                                        break;
                                case 5: neighbor[i][k]= base + (i+1)%l;//(right) right-up
                                        neighbor[i][k]=(neighbor[k][i]-l+sites)%sites;//(up)
                                        break;
                                case 6: neighbor[i][k]= base + (i+1)%l;//(right) right-down
                                        neighbor[i][k]= (neighbor[k][i]+l)%sites;//(down)
                                        break;
                                case 7: neighbor[i][k]= base + (i-1+l)%l;//(left)  left-down
                                        neighbor[i][k]= (neighbor[k][i]+l)%sites;//(down)
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
int randNeighbor(int id,int *vec,int id1,int id2,int size){
	int id_viz,idmax,idmin;

	if(id1<id2){
		idmin=id1;
		if(id2>=size)id2=size-1;
		idmax=id2;
	}else if(id2<id1){
		idmin=id2;
		if(id1>=size)id1=size-1;
		idmax=id1;
	}else{
		idmin=id1;
		idmax=idmin;
	}

	do{
		id_viz=(int)(FRANDOM*(idmax-idmin))+idmin;
	}while(id_viz==id);

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
int selectEvent(double *event,int size){
        int id;
        double nr,*cprob;

	cprob=(double *)calloc(size,sizeof(double));
	cumulProb(size,event,cprob);

        nr=FRANDOM*cprob[size-1];
        //id=bissectionSearch(nr, cprob, size);
        id=upperBoundStrict(nr, cprob, size);

	free(cprob);
        return id;
}
/***********************************************
*         randomly selects and event           *
*         using the respective cumulative      *
*         probability vector                   *
************************************************/
int selectEventCP(double *cprob,int size){
        int id;
        double nr;

        nr=FRANDOM*cprob[size-1];
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
 ***************************************************/
void logSpacedVec(double *vec,double vec0,double vecf,int size){
	int i;
	double x0,xi;

	x0=exp(vec0);
	vec[0]=vec0;
	vec[size-1]=vecf;
	for(i=1; i<size-1; ++i){
		xi=x0+(double)i*(exp(vecf)-exp(vec0))/size;
		vec[i]=log(xi);
	}

	return;
}
