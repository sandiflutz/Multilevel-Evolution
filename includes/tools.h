#ifndef TOOLS_H
#define TOOLS_H
/*************MACROS*****************************************************/
#define Pi 3.14159265358979323846  
#define max(a,b) (a>b?a:b)
#define min(a,b) (a>b?b:a)
/********structs*********************************************/
typedef struct{
        int *vec;//list of integers: vec[i]=id 
        int size;//real size
        int usize;//size of the part being used (usize<=size)
}DynList;//(fake) dynamic list

typedef struct{
	double *vecf;
	int sizef;
	int usizef;
}DynVec;
typedef struct{
	int sizeCL;
	int whichLB;
}ClusterMinimumID;
/********************************************************
 *                 Factorial                            *
 ********************************************************/
long int factorial(int n);
/********************************************************
*   Normal Distribution: calculates and returns         *
*   the probability density function value of a 
*   stardart normal in x				*
********************************************************/
double normalProb(double x);
/****************************************************************
*      Abramowitz and Stegun approximation to                   *
*      the complementary cumulative distribution function       *
*      for a standart normal (and Horner's method).             *
*      Source:                                                  *
*      https://www.johndcook.com/blog/normal_cdf_inverse/       *
*****************************************************************/
double rationalApprox(double t);
/****************************************************************
*     Return a value x from the inverse cumulative              *
*     distribution function of a standart normal                *
*****************************************************************/
double invNormalCDF(double p);
/****************************************************************
*     Return the cumulative distribution function for a normal  *
*     distribution of x                                         *
*****************************************************************/
double normalCDF(double x);
/****************************************************************
*     Generate a gaussian random number                         *
*****************************************************************/
double gaussRandNum(double mean, double sigma);
/*******************************************************************
*     Generate a gaussian random number (from a truncated dist.)   *
******************************************************************/
double truncGaussRandNum(double mean, double sigma,double xmin,double xmax);
/****************************************************************
*       Generate a sample of normaly distributed                *
*       random numbers in the range [xmin;xmax]                 *
*****************************************************************/
void normalRandSample(double mean, double sigma,double xmin,double xmax,int sample,double *randvec);
/****************************************************************
*      Draw a random integer from a poisson distribution o      *
*      mean @lambda                                             *
*****************************************************************/
int poissonRandKnuth(double lambda);
/************************************************************************************
*  Simple addition of an element to a list: add it at the end of the list and       *
*  increment the size of the list                                                   *
*************************************************************************************/
void listSimpleAdd(DynList *list,int add_elem);
/************************************************************************************
*  Simple subtraction of an element of a list: take the lest element of the         *
*  list and add it to the position of the element being subtracted and then         *
*  decrement the size of the list                                                   *
*************************************************************************************/
void listSimpleSub(DynList *list,int id_e);
/***********************************************************************************
*      Add an element of a list (size=@size and currently used size=@usize).       *
*      Elements added to the part of the list being used come from the             *
*      second part of the list (elements from @usize to @size-1)                   * 
*      Def: @add_elem is the element being add to the list,@id_e is the old        *
*      position of this element (in the second part of the list)                   *  
************************************************************************************/
void listAdd(DynList *list,int add_elem,int id_e);
/***********************************************************************************
*      Subtract an element of a list (size=@size and currently used size=@usize).  *
*      Elements subtracted from the list are stored in the second part of the list *
*      (elements from @usize to @size-1)                                           *
*      Def: @sub_elem is the element being subtracted from the list, @id_e is the  *
*      old position of this element (in the first part of the list)                *
************************************************************************************/
void listSub(DynList *list,int sub_elem,int id_e);
/****************************************************************
*       exchange 2 elements of a list of integers               *
*****************************************************************/
void exchange(int *list, int id1, int id2);
/****************************************************************
*       exchange 2 elements of a list of doubles                *
*****************************************************************/
void exchangeF(double *list, int id1, int id2);
/************************************************************************************
*    find the number of live neighbors of a specific host using the neighbors       *
*    matrix and the list of live hosts. Store live neighbors positions              *
************************************************************************************/
void searchLiveNeighbors(int sa,int idh,int *state,int **neighbors,DynList *alive_viz);
/************************************************************************************
*    find the number of empty sites in the neighborhood of a specific site using    *
*    the neighbors matrix. Store the neighboring empty sites positions              *
************************************************************************************/
void searchEmptyNeighbors(int se,int idh,int *state,int **neighbors,DynList *empty_viz);
/************************************************************************************
*    find the number of live neighbors using the neighbors matrix and              *
*    the list of live hosts. Store live neighbors positions                        *
*    If the nh*maxviz<nh/2 choose this instead (where nh is the number if live     * 
*    hosts ans maxviz in the number of neighbors)                                   *
************************************************************************************/
void searchAllLiveNeighbors(int dead,int sites,int nh,int maxviz,int *listh,int *host,int **neighbor,int *nviz,int *list_aviz);
/************************************************************************************
*    find the number of live neighbors using the links matrix and                   *
*    the list of live hosts. Store live neighbors positions.                        *
*    If the nh*maxviz>nh/2 choose this instead (where nh is the number if live      * 
*    hosts ans maxviz in the number of neighbors)                                   *
*************************************************************************************/
void searchAllLiveLinks(int sites,int nh,int *listh,int *netlink,int *nviz,int *list_aviz);
/*************************************************
* Set neighbors for a square lattice network     *
**************************************************/
void squareLattice(int **neighbor,int viz,int lsize);
/**********************************************************
*               Set network of a complete graph           *
*               (without self-neighboring)                *
***********************************************************/
void setCompleteGraph(int **neighbor,int sites);
/***********************************************************
*       Set a neighbors classification vector:             *
*       clneighbor[i][j]=k: j is the k-th neighbor of i    *
************************************************************/
void classifyNeighbors(int **neighbor,int **clneighbor,int nviz,int sites);
/********************************************************
*       Randomly select a neighbor from                 *
*       a neighbors list, from indexes 1d1 to id2-1     *
*********************************************************/
int randNeighborID(int id,int *vec,int id1,int id2,int size);
/*************************************************
*       Set network links between sites          *
**************************************************/
void setNetLinks(int **neighbor,int *netlink,int nviz,int sites);
/****************************************************
* Calculates the cumulative probability vector from *
* a simple probability vector                       *
*****************************************************/
void cumulProb(int size, double *prob,double *cprob);
/*********************************************************
 *   return the smallest index k that satisfy vec[k]>nr  *
 ********************************************************/
int upperBoundStrict(double nr,double *vec, int vsize);
/*********************************************************
*         randomly selects and event (towering method):  *
*   1)calculates the cumulative probabilitities          *
*   2)Randomly selects a number between 0 and the last   *
*   element of the cumulative probability vector         *
*   3)calls the bissectionSearch routine to find the     *
*   position, in the vector, the number belongs to.      *
*********************************************************/
int selectEvent(double nr,double *event,int size);
/***********************************************
*         randomly selects and event           *
*         using the respective cumulative      *
*         probability vector                   *
************************************************/
int selectEventCP(double nr,double *cprob,int size);
/***********************************************
*           Find the largest element           *
*           of a vector and returns            *
*           its position in the vector         *
************************************************/
int findMaxElement(int size,double *max_e,double *vec);
/********************************************************
 * return the smallest number  between @a and @b         *
 ********************************************************/
double minNum(double a, double b);
/********************************************************
 *	Create fill a vector with log spaced elements	*
 ********************************************************/
void logSpacedVec(double *vec,double e0,double ef,double m0, double mf, int size);
/********************************************************
 * 		Build a histogram			*
 ********************************************************/
void buildHistogram(double binsize,int nbins,int visize,double *vec,int *hist);
/********************************************************
*  calculate the root mean square of the elements       *
*  from index idi to index idf of a vector vec[]        *
*********************************************************/
void calcRMSError(double *vec,int idi, int idf,double *stats);
/***************************************************
*  calculate spatial corretation for a square      *
*  lattice type of system                          *
****************************************************/
double spatialCorr(int *state,int sites, int dist,int right,int down,int **neighbor);
/***************************************************
*  calculate spatial corretation for a square      *
*  lattice type of system in 1 direction only      *
****************************************************/
double spatialCorr1d(int *state,int sites, int dist,int id_direction,int **neighbor);
/****************************************
 *      Heaviside Step Function         *
 *      for x in (-1.,1]                *
 ***************************************/
double stepFuncBounded(double x);
/****************************************************************
*       set cluster labels (Hoshen-Kopelman algorithm)          *
*****************************************************************/
int setClLabelsHK(int *state,int **neighbor,int vleft,int vup,int *labels,int lsize);
/****************************************************************
*       set cluster labels (Hoshen-Kopelman algorithm)          *
*       using a list of occupied sites and the inverse list     *
*       -list of occupied sites, list, is a DynList struct      *
*       (that contains a vector of integers,vec[], an integer   *
*       called usize=the number of occupied sites, and an       *
*       integer called size=the number of sites)                *
*       -the first list->usize positions of the list->vec[]     *
*       receive the positions of the occupied sites             *
*       -the last (list->size)-(list->usize) positions receive  *
*       the positions of the unoccupied sites                   *
*       -the vector ilist is the inverse list: its indexes are  *
*       the site positions, that go from 0 to list->size, and   *
*       its elements are the indexes of those positions on the  *
*       list of occupied sites, list->vec                       *
*****************************************************************/
void setClusterLabelsWithList(DynList *list,int *ilist,int **neighbor,int vleft,int vup,int *labels);
/********************************************************
*       Part of the Hoshen-Kopelman algorithm.          *
*********************************************************/
int find(int x,int *lblist);
/********************************************************
*	Part of the Hoshen-Kopelman algorithm.		*
*	Update the list of the list of labels, lblist:	*
*	the position with the largest label between 	*
*	x and y, in @lblist, stores the smalest label.	*
*	Return the smallest label.			*
*********************************************************/
int unionFind(int x,int y,int *lblist);
/****************************************************************
*  Fix cluster label order so labels are =0,1,...,n-1,          *
*  with n=number of clusters. To be used, if necessary, after   *
*  the Hoshen-Kopelman algorithm (with the least label of a     *
*  cluster being equal to the least node id in the cluster).    *
*  Returns the number of clusters.                              *
*****************************************************************/
int fixClusterLbOrder(int nid,int *labels);
/************************************************************************
*       Measure the average cluster size, the related standart          *
*       deviation and the size and label of the largest cluster.        *
*************************************************************************/
void calcClusterSizeStats(int nid,int ncl,int *labels,int *clsize,double *stats,ClusterMinimumID *maxclid);
#endif
