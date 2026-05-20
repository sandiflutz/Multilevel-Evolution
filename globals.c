/*globals.c*/
#include"tools.h"
#include"globals.h"
int *host = NULL;
int **neighbor = NULL;
int *inverselisth = NULL;
int *lb = NULL;
ClusterFullID *maxclw =NULL;
double *rho_e = NULL;
double *bac = NULL;
double **costvec;
double *dtVec = NULL;
DynList *listh = NULL;
DynVec *offcomp = NULL;
DynVec *avinv = NULL;
SysParams *spar = NULL;
SysTimes *stime = NULL;
GenFile *gfile = NULL;
EvMeasures *meas = NULL;
SysMeasures *sysmeas = NULL;
DynList *bestwlisth = NULL;
