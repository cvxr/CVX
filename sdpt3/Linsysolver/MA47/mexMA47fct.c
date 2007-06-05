/*************************************************************************/
/* Usage: [A, IW, Flag, Info] = ma47fct(A, Keep, JCN, La, Liw, Cntl, ICtnl); */
/*        [A, IW, Flag, Info] = ma47fct(A, Keep, JCN, La, Liw, Cntl); */
/*        [A, IW, Flag, Info] = ma47fct(A, Keep, JCN, La, Liw); */
/*************************************************************************/

#if !defined(max)
#define  max(A, B)   ((A) > (B) ? (A) : (B))
#endif

#include "mex.h"
#include "MA47.h"

void convertin(int n, int lkeep, int ne, int *mxir, int *mxjc,int *keep, int *jcn, 
	       double *mxA, double *A, double *mxkeep, double *mxjcn);
void cpoutput(int la, int liw, int linfo, int *iw, int *info, double *A, double *mxA,
	      double *mxiw, double *mxinfo);
int checkinfo(int *info, int *tmpinfo, int *flag);

/*************************************************************************/
void mexFunction(const int nlhs, mxArray *plhs[],
		 const int nrhs, const mxArray *prhs[])
{
    const mxArray  *Lsymb; 
    int    i;
    int    n, ne, liw, lkeep, la, linfo, liw1;
    int    *piw, *pinfo,  *picntl, *pkeep, *pjcn, *pflag;
    int    *pmxjc, *pmxir, *piw1, *tmpinfo;
    double *pmxicntl, *pmxcntl,*pmxiw, *pmxinfo, *pmxA, *pmxjcn, *pmxflag, *pmxkeep; 
    double *prinfo, *pcntl, *pA, *pnz;

/* KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS */
/* FOR USE IN ALL YOUR FORTRAN MEX FILES. */
/* ---------------------------------------------------------------*/


/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

    /* Function Body */
    if (nrhs < 2) {
       mexErrMsgTxt("mexMA47fct requires at least 2 input arguments");
    } else if (nlhs < 2) {
       mexErrMsgTxt("mexMA47fct requires at least 2 output arguments");
    }

/* Verify input size  */

    n = mxGetM(prhs[0]);
    if (!mxIsSparse(prhs[0])){
       mexErrMsgTxt("mexMA47fct:  A must be sparse"); }
    if (n != mxGetN(prhs[0])) {
       mexErrMsgTxt("mexMA47fct:  A must be square"); }
    pmxir   = mxGetIr(prhs[0]);
    pmxjc   = mxGetJc(prhs[0]);
    pmxA    = mxGetPr(prhs[0]);

/* SPECIFY THE DIMENSIION OF WORKING VECTORS */

    if(!mxIsStruct(prhs[1])) { 
       mexErrMsgTxt("mexMA47fct: 2nd input should be a structure."); }    
    if( (Lsymb = mxGetField(prhs[1],0,"keep")) == NULL) {
       mexErrMsgTxt("Missing field Lsymb.keep"); }
    lkeep = max(mxGetM(Lsymb),mxGetN(Lsymb));
    pmxkeep = mxGetPr(Lsymb);
    if( (Lsymb = mxGetField(prhs[1],0,"Jcn")) == NULL) {
       mexErrMsgTxt("Missing field Lsymb.Jcn"); }
    ne = max(mxGetM(Lsymb),mxGetN(Lsymb));
    pmxjcn = mxGetPr(Lsymb);
    if( (Lsymb = mxGetField(prhs[1],0,"La")) == NULL) {
       mexErrMsgTxt("Missing field Lsymb.La"); }    
    la    = (int) mxGetScalar(Lsymb);
    if( (Lsymb = mxGetField(prhs[1],0,"Liw")) == NULL) {
       mexErrMsgTxt("Missing field Lsymb.Liw"); }    
    liw   = (int) mxGetScalar(Lsymb);
    liw1  = 2*n+2;

/* Set paramenters in CNTL & ICNTL  */

    picntl = (int *)mxCalloc(7,    sizeof(int));
    pcntl  = (double *)mxCalloc(2,  sizeof(double));
    ma47id_(pcntl, picntl);
    pcntl[0] = 0.00001;

    if (nrhs > 2){
       if (max(mxGetM(prhs[2]),mxGetN(prhs[2])) != 2){
	  mexErrMsgTxt("mexMA47fct: CNTL requires 2 parameters"); }
       pmxcntl = mxGetPr(prhs[2]);
       pcntl[0] = pmxcntl[0];
       pcntl[1] = pmxcntl[1];
    }
    if (nrhs > 3){
       if (max(mxGetM(prhs[3]),mxGetN(prhs[3])) != 7){
	   mexErrMsgTxt("mexMA47fct: ICNTL requires 7 parameters"); }
       pmxicntl = mxGetPr(prhs[3]);
     /* convert the parameter form Matlab double format to integer*/
       for (i=0; i<7; i++) { picntl[i] = (int)pmxicntl[i]; }
    }

/* CREATE WORKING PARAMETERS */

    pjcn    = (int *)mxCalloc(ne,    sizeof(int));
    pkeep   = (int *)mxCalloc(lkeep, sizeof(int));
    piw     = (int *)mxCalloc(liw,   sizeof(int));
    pinfo   = (int *)mxCalloc(24,    sizeof(int));
    piw1    = (int *)mxCalloc(liw1,  sizeof(int));
    prinfo  = (double *)mxCalloc(4,  sizeof(double));
    pnz     = (double *)mxCalloc(ne, sizeof(double));

/* INPUT DATA TRANSFORMATION */
    convertin(n,lkeep,ne,pmxir,pmxjc,pkeep,pjcn,pmxA,pnz,pmxkeep,pmxjcn);

/* DO THE ACTUAL COMPUTATIONS IN A FORTRAN SUBROUTINE */
    pA  = (double *)mxCalloc(la, sizeof(double));
    for(i=0; i<ne; i++) { pA[i] = pnz[i]; }
    ma47bd_(&n,&ne,pjcn,pA,&la,piw,&liw,pkeep,pcntl,picntl,piw1,prinfo,pinfo);

/* CHECK ERROR FLAG */
      while(pinfo[0] == -3 || pinfo[0] == -4){
	   if (pinfo[0] == -3) {
	      liw  = (int)(1.2*pinfo[1]);
	      piw  = (int *)mxRealloc(piw, liw*sizeof(int));
	      ma47bd_(&n,&ne,pjcn,pA,&la,piw,&liw,pkeep,pcntl,picntl,piw1,prinfo,pinfo);
	   }
	   if (pinfo[0] == -4) {
	      la = (int)(1.2*pinfo[1]);
	      pA = (double *)mxRealloc(pA, la*sizeof(double));
	      for(i=0; i<ne; i++){
                 pA[i] = pnz[i]; }
	      ma47bd_(&n,&ne,pjcn,pA,&la,piw,&liw,pkeep,pcntl,picntl,piw1,prinfo,pinfo);
	   }
      }
      tmpinfo = (int *)mxCalloc(1,  sizeof(int));
      pflag   = (int *)mxCalloc(1,  sizeof(int));
      linfo   = checkinfo(pinfo, tmpinfo, pflag);
	  
/* CREATE OUTPUT PARAMETERS */

      plhs[0] = mxCreateDoubleMatrix(la, 1, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(liw, 1, mxREAL);
      plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
      plhs[3] = mxCreateDoubleMatrix(linfo, 1, mxREAL);

/* DEREFERENCE OUTPUT ARGUMENTS TO GET REAL PART POINTERS */

      pmxA    = mxGetPr(plhs[0]);
      pmxiw   = mxGetPr(plhs[1]);
      pmxflag = mxGetPr(plhs[2]);
      pmxinfo = mxGetPr(plhs[3]);      
     *pmxflag = *pflag;         
      cpoutput(la,liw,linfo,piw,tmpinfo,pA,pmxA,pmxiw,pmxinfo);

/* Release working arrays */

	mxFree(picntl);
	mxFree(pcntl);
	mxFree(pjcn);
	mxFree(pnz);
	mxFree(pkeep);
	mxFree(piw);
	mxFree(prinfo);
	mxFree(pinfo);
	mxFree(piw1);
	mxFree(tmpinfo);
	mxFree(pA);
	mxFree(pflag);

} /* mexFunction*/


/* ---------------------------------------------------------------- */
/* Convert from Matlab sparse format to HSL sparse format  */
/* ---------------------------------------------------------------- */

void convertin(int n, int lkeep, int ne, int *mxir, int *mxjc,int *keep, int *jcn, 
	       double *mxA, double *A, double *mxkeep, double *mxjcn)

{   int i, j, k;

	for (i=0, k=0; i<n; i++) {
	    for (j = mxjc[i]; j< mxjc[i+1]; j++){
		if (mxir[j] >= i) { A[k++] = mxA[j]; }		  
	    }
	}
	for (i=0; i<lkeep; ++i) {
	    keep[i] = (int)mxkeep[i];
	}
	for (i=0; i<ne; ++i) {
	    jcn[i] = (int)mxjcn[i];
	}

} /* convertin*/

/* ------------------------------------------------- */
/* Convert outputs from integer type in Fortran */
/*                   to real*8  type in Matlab */
/*------------------------------------------------- */
void cpoutput(int la, int liw, int linfo, int *iw, int *info, double *A, double *mxA,
	      double *mxiw, double *mxinfo)

{   int i;

    for (i = 0; i<la; i++) {
	mxA[i] = A[i];
    }
    for (i = 0; i < liw; ++i) {
	mxiw[i] = (double) iw[i];
    }
    for (i = 0; i < linfo; ++i) {
	mxinfo[i] = (double) info[i];
    }

}/* cpoutput */

/* ----------------------------------------------------------------- */
/* Check the infomation given by ma47bd_(...). */
/* In the case of error, extract the corresponding infomation for output */
/* ----------------------------------------------------------------- */
int checkinfo(int *info, int *tmpinfo, int *flag)

{  int n;

   *flag = info[0];
    n = 0;
    if (info[0] == 4){
       n = 1;
       tmpinfo[0] = info[23];
    }
    return n;
}
/* ----------------------------------------------------------------- */
