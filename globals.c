/*globals.c*/
#include"tools.h"
#include"globals.h"
int *host = NULL;
int **neighbor = NULL;
int *netlink = NULL;
double **bac = NULL;
double *dtVec = NULL;
char *fdatapath;
DynList *listh = NULL;
DynList *listnb = NULL;
DynList *alive_viz = NULL;
SysParams *spar = NULL;
