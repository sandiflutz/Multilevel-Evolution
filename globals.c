/*globals.c*/
#include"tools.h"
#include"globals.h"
int *host = NULL;
int **neighbor = NULL;
int **clneighbor = NULL;
int *netlink = NULL;
int *inverselisth=NULL;
double **bac = NULL;
double *dtVec = NULL;
char *fdatapath;
DynList *listh = NULL;
Dyn2DList *listh2d = NULL;
DynList *alive_viz = NULL;
SysParams *spar = NULL;
