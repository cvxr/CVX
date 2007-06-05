/*************************************************************************/
/* Usage: [Sol] = ma47slv(A, IW, Rhs, ICntl); */
/*        [Sol] = ma47slv(A, IW, Rhs); */ 
/*************************************************************************/

#if !defined(max)
#define  max(A, B)   ((A) > (B) ? (A) : (B))
#endif

#include "mex.h"
#include "MA47.h"

/* ---------------------------------------------------------------------*/
/* KEEP THE FOLLOWING SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS */
/* FOR USE IN ALL YOUR FORTRAN MEX FILES. */
/* ---------------------------------------------------------------------*/

void convertin(int liw, double *mxiw, int *iw);
void cpoutput(int n, double *mxrhs, double *mxsol);

/*************************************************************************/

void mexFunction(const int nlhs, mxArray *plhs[],
		 const int nrhs, const mxArray *prhs[])
{	
    int    n, la, liw, i, issprhs, nzrhs, k;
    int    *piw, *picntl, *piw1, *rhsIr, *rhsJc; 
    double *pmxsol, *pmxiw, *pmxicntl, *pmxA, *pmxrhs, *rhstmp; 
    double *pw, *pcntl;

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

    if (nrhs < 3) {
       mexErrMsgTxt("mexMA47slv requires at least 3 input arguments");
    } else if (nlhs != 1) {
       mexErrMsgTxt("mexMA47slv requires 1 output arguments");
    }

    /* Set paramenters in ICNTL */
    picntl = (int *)mxCalloc(7,    sizeof(int));
    pcntl  = (double *)mxCalloc(2, sizeof(double));
    ma47id_(pcntl, picntl);

    if (nrhs >3){
       if (max(mxGetM(prhs[3]),mxGetN(prhs[3])) != 7){
	  mexErrMsgTxt("mexMA47slv: ICNTL requires 7 parameters");
       }
       pmxicntl = mxGetPr(prhs[3]);
       /* convert the parameter form Matlab double format to integer*/
       for (i=0; i<7; i++){
	   picntl[i] = (int)pmxicntl[i]; }		
    }

/* SPECIFY THE DIMENSIION OF WORKING VECTORS */	

    la  = max(mxGetM(prhs[0]),mxGetN(prhs[0]));
    liw = max(mxGetM(prhs[1]),mxGetN(prhs[1]));
    n   = max(mxGetM(prhs[2]),mxGetN(prhs[2])); 

/* DEREFERENCE ARGUMENTS TO GET ARRAY POINTERS */

    pmxA   = mxGetPr(prhs[0]);
    pmxiw  = mxGetPr(prhs[1]);
    pmxrhs = mxGetPr(prhs[2]); 
    issprhs = mxIsSparse(prhs[2]);
    rhstmp  = mxCalloc(n,sizeof(double));   
    if (!issprhs) {
       for (k=0; k<n; k++) {
	   rhstmp[k] = pmxrhs[k]; }
    }
    else {
       rhsJc = mxGetJc(prhs[2]); 
       rhsIr = mxGetIr(prhs[2]); 
       nzrhs = rhsJc[1]; 
       for (k=0; k<nzrhs; k++) {
	   rhstmp[rhsIr[k]] = pmxrhs[k]; }
    }
/* CREATE WORKING PARAMETERS */

    piw  = (int *)mxCalloc(liw,  sizeof(int));
    piw1 = (int *)mxCalloc(n,    sizeof(int));
    pw   = (double *)mxCalloc(n, sizeof(double));

/* INPUT DATA TRANSFORMATION */

    convertin(liw,pmxiw,piw);

/* DO THE ACTUAL COMPUTATIONS IN A FORTRAN SUBROUTINE */

    ma47cd_(&n,pmxA,&la,piw,&liw,pw,rhstmp,piw1,picntl); 

/* CREATE OUTPUT PARAMETERS */

    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    pmxsol  = mxGetPr(plhs[0]);

/* OUTPUT DATA TRANSFORMATION */

    cpoutput(n,rhstmp,pmxsol);
      
/* Release working arrays */

    mxFree(rhstmp); 
    mxFree(picntl);
    mxFree(pcntl);
    mxFree(piw);
    mxFree(pw);
    mxFree(piw1);

} /* mexFunction*/

/* ---------------------------------------------------------------- */
/* Convert from type real*8 in Matlab to type integer in Fortran */
/* ---------------------------------------------------------------- */
void convertin(int liw, double *mxiw, int *iw)

{   int i;

    for (i=0; i < liw; ++i) {
        iw[i] = (int) mxiw[i];
    }
} /* convertin*/

/* ---------------------------------------------------------------- */
/* Copy from prhs[*]* to plhs[*]*        */
/* ---------------------------------------------------------------- */
void cpoutput(int n, double *mxrhs, double *mxsol)

{   int i;

    for (i = 0; i < n; ++i) {
	mxsol[i] = mxrhs[i];
    }
} /* cpoutput*/
/* ---------------------------------------------------------------- */
