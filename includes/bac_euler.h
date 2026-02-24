#ifndef BAC_EULER_H
#define BAC_EULER_H

/********************************************
 * integration of the microbial equations   *
 * using euler method                       *
 ********************************************/
void bac_euler(double dt,SysParams *sp);
/****************************************************************
 *     cost reduction function caused by negative types         *
 *     of bacteria                                              *
 ***************************************************************/
double expCostReduc(double x, double a,double b,double eps);
#endif

