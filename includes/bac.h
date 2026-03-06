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
/********************************************
 * integration of the microbial equations   *
 * using euler method                       *
 ********************************************/
void bacDyn(int idh,double dt,DynList *liveviz,double **btmp);
/****************************************************************
 * Microbial layer evolution (for a time interval=dt)   *
 ****************************************************************/
void evoBac(double dt);
#endif

