/*globals.c*/
#include"tools.h"
#include"globals.h"
int *host = NULL;
int **neighbor = NULL;
int **clneighbor = NULL;
int *netlink = NULL;
int *inverselisth=NULL;
double **bac = NULL;
double **costvec;
double *dtVec = NULL;
DynList *listh = NULL;
DynList *alive_viz = NULL;
DynVec *offcomp = NULL;
DynVec *avinv = NULL;
SysParams *spar = NULL;
SysTimes *stime = NULL;
GenFile *gfile = NULL;
EvMeasures *meas = NULL;
