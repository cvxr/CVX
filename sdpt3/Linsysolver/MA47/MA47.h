

#ifndef MA47_H
#define MA47_H

int ma47ad_(int *n, int *ne, int *irn, int *jcn, int *iw, int *liw, 
	    int *keep, int *icntl, double *rinfo, int *info);
int ma47bd_(int *n, int *ne, int *jcn, double *a, int *la, int *iw, 
	    int *liw, int *keep, double *cntl, int *icntl, int *iw1, 
	    double *rinfo, int *info);
int ma47cd_(int *n, double *a, int *la, int *iw, int *liw, double *w, 
	    double *rhs, int *iw1, int *icntl);
int ma47id_(double *cntl, int *icntl);


#endif

