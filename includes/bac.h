#ifndef BAC_H
#define BAC_H
/****************************************************************
 *     cost reduction function caused by negative types         *
 *     of bacteria                                              *
 ***************************************************************/
double expCostReduc(double x, double a,double b,double eps);
/************************************************************************
*      update costs multiplicative factor acording to the amount       *
*      of negative microbial types                                     *
************************************************************************/
void updateCosts(int idh,double mtot,double c0,double *w,double **cvec,double **bstate);
/********************************************************
*       integration of the microbial equations          *
*       using euler method when the cost of helping     *
*       depends bacteria abundances                     *
 ********************************************************/
void bacDynCostVec(int idh,double dt,double *btmp,double time,double *type_abund,int nv);
/********************************************
 * integration of the microbial equations   *
 * using euler method                       *
 ********************************************/
void bacDyn(int idh,double dt,double *btmp,double time,double *type_abund,int nv);
/****************************************************************
 * Microbial layer evolution (for a time interval=dt)   *
 ****************************************************************/
void evoBac(double dt,double time);
#endif

