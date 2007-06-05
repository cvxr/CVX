/***********************************************************************
* mexnnz.c: get number of non-zero elements of a matrix M 
*
*  nnz = mexnnzt(M); 
*
* SDPT3: version 3.0
* Copyright (c) 1997 by
* K.C. Toh, M.J. Todd, R.H. Tutuncu
* Last Modified: 2 Feb 01
***********************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

/**********************************************************
* 
***********************************************************/
void mexFunction(int nlhs, mxArray  *plhs[], 
                 int nrhs, const mxArray  *prhs[] )

{    double   *A, *nnz;  
     mwIndex  *irA, *jcA; 

     int      NZmax, m, n, isspA, j, k, jm, nnztmp;

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

     if (nrhs != 1){
         mexErrMsgTxt("mexnnz: requires 1 input arguments."); }
     else if (nlhs>1){ 
         mexErrMsgTxt("mexnnz: requires 1 output argument."); }

/* CHECK THE DIMENSIONS */

     plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
     nnz = mxGetPr(plhs[0]); 
     if (mxGetPr(prhs[0]) == NULL) { 
        nnz[0] = 0.0; 
        return; 
     }   
     m = mxGetM(prhs[0]); 
     n = mxGetN(prhs[0]);
     A = mxGetPr(prhs[0]);
     isspA = mxIsSparse(prhs[0]); 
     if (isspA) { 
        irA = mxGetIr(prhs[0]);      
        jcA = mxGetJc(prhs[0]);
        NZmax = jcA[n]; 
     }
/***** main body *****/ 
     nnztmp = 0; 
     if (isspA) { 
        for (k=0; k<NZmax; k++) { 
	   if (A[k] != 0) { nnztmp++; }
        }
     } else {
        for (j=0; j<n; j++) { 
           jm = j*m; 
	   for (k=0; k<m; k++) {
	      if (A[k+jm] !=0) { nnztmp++; } 
	   }
        }
     }      
     nnz[0] = (double)nnztmp;
 return;
 }
/**********************************************************/

