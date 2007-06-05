/*************************************************************************/
/* Usage: [Keep, JCN, Flag, Info] = ma47syb(NE, A, ICntl) */
/*        [Keep, JCN, Flag, Info] = ma47syb(NE, A) */
/*************************************************************************/

#if !defined(max)
#define  max(A, B)   ((A) > (B) ? (A) : (B))
#endif

#include "mex.h"
#include "MA47.h"

void convertin(int n, int *mxirn, int *mxjcn, int *irn, int *jcn);
void cpoutput(int lkeep, int ne, int linfo, int *keep, int *jcn, int *info, 
	      double *mxkeep, double *mxjcn, double *mxinfo);
int checkinfo(int *info, int *tmpinfo, int *flag);

void mexFunction(const int nlhs, mxArray *plhs[],
		 const int nrhs, const mxArray *prhs[])
{
    /* Local variables */
    int n, ne, liw, lkeep,linfo, i;
    int *piw, *pinfo,  *picntl, *pkeep, *pjcn,  *pirn, *pflag;
    int *pmxjcn, *pmxirn, *tmpinfo;
    double *pmxicntl, *pmxkeep, *pmxinfo, *pmxflag, *pmxjcn_d; 
    double *prinfo, *pcntl;

/* KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS */
/* FOR USE IN ALL YOUR FORTRAN MEX FILES. */
/* --------------------------------------------------------------------- 
*/


/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

    /* Function Body */
    if (nrhs < 2) {
       mexErrMsgTxt("mexMA47syb requires at least 2 input arguments");
    } else if (nlhs != 4) {
       mexErrMsgTxt("mexMA47syb requires 4 output arguments");
    }
   
    /* Set paramenters in ICNTL */
    picntl = (int *)mxCalloc(7,    sizeof(int));
    pcntl  = (double *)mxCalloc(2,  sizeof(double));
    ma47id_(pcntl, picntl);

    if (nrhs > 2){
       if (max(mxGetM(prhs[2]),mxGetN(prhs[2])) != 7){
	  mexErrMsgTxt("mexMA47syb: ICNTL requires 7 parameters");
       }
       pmxicntl = mxGetPr(prhs[2]);
       /* convert the parameter form Matlab double format to integer*/
       for (i=0; i<7; i++) { picntl[i] = (int)pmxicntl[i]; }	 
    }

/* verify 2nd input */
    n = mxGetM(prhs[1]);
    if (!mxIsSparse(prhs[1])){
	mexErrMsgTxt("mexMA47syb:  A must be sparse"); }
    if (n != mxGetN(prhs[1])) {
	mexErrMsgTxt("mexMA47syb:  A must be square"); }

/* SPECIFY THE DIMENSIION OF WORKING VECTORS */
    ne     = (int) mxGetScalar(prhs[0]);;
    liw    = (int)(((ne << 1) + n*5 + 4) * 1.2);
    lkeep  = ne + 5*n + 2;

/* DEREFERENCE ARGUMENTS TO GET ARRAY POINTERS */
    pmxirn  = mxGetIr(prhs[1]);
    pmxjcn  = mxGetJc(prhs[1]);

/* CREATE WORKING PARAMETERS */
    pirn    = (int *)mxCalloc(ne,    sizeof(int));
    pjcn    = (int *)mxCalloc(ne,    sizeof(int));
    pkeep   = (int *)mxCalloc(lkeep, sizeof(int));
    piw     = (int *)mxCalloc(liw,   sizeof(int));
    pinfo   = (int *)mxCalloc(24,    sizeof(int));
    prinfo  = (double *)mxCalloc(4,  sizeof(double));

/* INPUT DATA TRANSFORMATION */

    convertin(n, pmxirn, pmxjcn, pirn, pjcn);

/* DO THE ACTUAL COMPUTATIONS IN A FORTRAN SUBROUTINE */

    ma47ad_(&n,&ne,pirn,pjcn,piw,&liw,pkeep,picntl,prinfo,pinfo);

/* CHECK ERROR FLAG */

    while (pinfo[0] == -3) {
        liw  = (int)(1.2*pinfo[1]);
        piw  = (int *)mxRealloc(piw, liw*sizeof(int));
        ma47ad_(&n,&ne,pirn,pjcn,piw,&liw,pkeep,picntl,prinfo,pinfo);
    }
    tmpinfo = (int *)mxCalloc(3,  sizeof(int));
    pflag = (int *)mxCalloc(1,  sizeof(int));
    linfo = checkinfo(pinfo, tmpinfo, pflag);

/* CREATE OUTPUT PARAMETERS */

      plhs[0]  = mxCreateDoubleMatrix(lkeep, 1, mxREAL);
      plhs[1]  = mxCreateDoubleMatrix(ne, 1, mxREAL);
      plhs[2]  = mxCreateDoubleMatrix(1, 1, mxREAL);
      plhs[3]  = mxCreateDoubleMatrix(linfo, 1, mxREAL);

/* DEREFERENCE OUTPUT ARGUMENTS TO GET REAL PART POINTERS */

      pmxkeep  = mxGetPr(plhs[0]);     
      pmxjcn_d = mxGetPr(plhs[1]); 
      pmxflag  = mxGetPr(plhs[2]);
      pmxinfo  = mxGetPr(plhs[3]);
     *pmxflag  = (double)*pflag;

      cpoutput(lkeep,ne,linfo,pkeep,pjcn,tmpinfo,pmxkeep,pmxjcn_d,pmxinfo);

/* Release working arrays */

	mxFree(picntl);
	mxFree(pcntl);
	mxFree(pirn);
	mxFree(pjcn);
	mxFree(pkeep);
	mxFree(piw);
	mxFree(prinfo);
	mxFree(pinfo);
	mxFree(tmpinfo);
	mxFree(pflag);


} /* mexFunction*/

/* ---------------------------------------------------------------- */
/* Convert from Matlab sparse format to HSL sparse format  */
/* ---------------------------------------------------------------- */

void convertin(int n, int *mxirn, int *mxjcn, int *irn, int *jcn)

{  int i, j, k;

    for (i = 0, k = 0; i < n; i++) {
	for (j = mxjcn[i]; j< mxjcn[i+1]; j++){
	    if (mxirn[j] >= i){
	       jcn[k] = i+1;
	       irn[k] = mxirn[j]+1;
	       k++;
	    }
	}
    }
} /* convertin*/


/* ------------------------------------------------- */
/* Convert outputs from integer type in Fortran */
/*                   to real*8  type in Matlab */
/*------------------------------------------------- */
void cpoutput(int lkeep, int ne, int linfo, int *keep, int *jcn, int *info, 
	      double *mxkeep, double *mxjcn, double *mxinfo)

{   int i;

    for (i = 0; i < lkeep; ++i) {
	mxkeep[i] = (double) keep[i];
    }
     for (i = 0; i < ne; ++i) {
	 mxjcn[i] = (double) jcn[i];
    }
    for (i = 0; i < linfo; ++i) {
	mxinfo[i] = (double) info[i];
    }	
}/* cpoutput */


/* ----------------------------------------------------------------- */
/* Check the infomation given by ma47ad_(...). */
/* In the case of error, extract the corresponding infomation for output */
/* ----------------------------------------------------------------- */
int checkinfo(int *info, int *tmpinfo, int *flag)

{	int n;

	*flag = info[0];

	switch(info[0]){
	case -5:
	case -6:
	    n = 1;
	    tmpinfo[0] = info[1];
	    break;
	case 0:
 	    n = 2;
	    tmpinfo[0] = info[5];
	    tmpinfo[1] = info[6];
	    break;
	case 1:
	    n = 1;
	    tmpinfo[0] = info[2];
	    break;
	case 2:
	    *flag = 0;
	    n = 3;
	    tmpinfo[0] = info[5];
	    tmpinfo[1] = info[6];
	    tmpinfo[2] = info[3];
	    break;
	case 3:
	    n = 2;
	    tmpinfo[0] = info[2];
	    tmpinfo[1] = info[3];
	    break;
	case 4:
	    n = 1;
	    tmpinfo[0] = info[7];
	    break;
	case 5:
	    n = 2;
	    tmpinfo[0] = info[2];
	    tmpinfo[1] = info[7];
	    break;
	case 6:
	    n = 2;
	    tmpinfo[0] = info[3];
	    tmpinfo[1] = info[7];
	    break;
	case 7:
	    n = 3;
	    tmpinfo[0] = info[2];
	    tmpinfo[1] = info[3];
	    tmpinfo[2] = info[7];
	    break;
	default:
	    n = 0;
	    break;
	}
	return n;
}
/* ----------------------------------------------------------------- */
