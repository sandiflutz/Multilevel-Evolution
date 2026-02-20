#ifndef TOOLS_H
#define TOOLS_H
/*************MACROS*****************************************************/
#define Pi 3.14159265358979323846  
/********structs*********************************************/
typedef struct{
        int *vec;//list of integers: vec[i]=id 
        int size;//real size
        int usize;//size of the part being used (usize<=size)
}DynList;//(fake) dynamic list

typedef struct{
        double *vecf;//list of integers
        int sizef;//real size
        int usizef;//size of the part being used (usize<=size)
}DynListF;//(fake) dynamic list of fluctuating point numbers 

typedef struct{
        double **array;//list of integers
        int ncol;//# number of columns
        int nlin;//# number of lines
        int *ucol;//# of columns being used 
        int ulin;//# of lines being used 
}Dyn2DList;//(fake) 2D dynamic list of fluctuating
/**************************************************************/
/********************************************************
 *                 Factorial                            *
 ********************************************************/
long int factorial(int n);
/********************************************************
*   Normal Distribution: calculates and returns         *
*   the probability of x in a normal distribution       *
*   with standard deviation of 1 and mean of 0          *
********************************************************/
double normalProb(double x);
/****************************************************************
*     Generate a gaussian random number                         *
*****************************************************************/
double gaussRandNum(double mean, double var);
/*******************************************************************
*     Generate a gaussian random number (from a truncated dist.)   *
******************************************************************/
double truncGaussRandNum(double mean, double var,double a,double b);
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
/**************Double List Version***********************************/
void listSimpleAddF(DynListF *list,double add_elem);
/************************************************************************************
*  Simple subtraction of an element of a list: take the lest element of the         *
*  list and add it to the position of the element being subtracted and then         *
*  decrement the size of the list                                                   *
*************************************************************************************/
void listSimpleSub(DynList *list,int id_e);
/**************Double List Version***********************************/
void listSimpleSubF(DynListF *list,int id_e);
/***********************************************************************************
*      Add an element of a list (size=@size and currently used size=@usize).       *
*      Elements added to the part of the list being used come from the             *
*      second part of the list (elements from @usize to @size-1)                   * 
*      Def: @add_elem is the element being add to the list,@id_e is the old        *
*      position of this element (in the second part of the list)                   *  
************************************************************************************/
void listAdd(DynList *list,int add_elem,int id_e);
/**************Double List Version***********************************/
void listAddF(DynListF *list,double add_elem,int id_e);
/***********************************************************************************
*      Subtract an element of a list (size=@size and currently used size=@usize).  *
*      Elements subtracted from the list are stored in the second part of the list *
*      (elements from @usize to @size-1)                                           *
*      Def: @sub_elem is the element being subtracted from the list, @id_e is the  *
*      old position of this element (in the first part of the list)                *
************************************************************************************/
void listSub(DynList *list,int sub_elem,int id_e);
/**************Double List Version***********************************/
void listSubF(DynListF *list,double sub_elem,int id_e);
/************************************************************
*               exchange 2 elements of a list               *
*************************************************************/
void exchange(int *list, int id1, int id2);
/************************************************************************************
*    find the number of live neighbors of a specific host using the neighbors       *
*    matrix and the list of live hosts. Store live neighbors positions              *
************************************************************************************/
void searchLiveNeighbors(int dead,int idh,int *host,int **neighbor,DynList *alive_viz);
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
void squareLattice(int **neighbor,int viz,int sites);
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
/**************************************************
*       Randomly select a neighbor from           *
*       a neighbors list, from indexes 1d1 to id2 *
***************************************************/
int randNeighbor(int id,int *vec,int id1,int id2,int size);
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
int selectEvent(double *event,int size);
/***********************************************
*         randomly selects and event           *
*         using the respective cumulative      *
*         probability vector                   *
************************************************/
int selectEventCP(double *cprob,int size);
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
/****************************************************
 * Create fill a vector with log spaced elements    *
 ***************************************************/
void logSpacedVec(double *vec,double e0,double ef,double m0, double mf, int size);
#endif
