/***********************************************************************
* mexexpand.c : C mex file 
*
*   z = mexexpand(blk,x); 
* 
*   Input: blk   = [n1, n2, ..., nk]
*
*   Output: [x1 x1...x1,  x2 x2...x2, ...., xk xk...xk]'    
*            n1          n2                nk
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

#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

/**********************************************************
* 
***********************************************************/
void mexFunction(int nlhs, mxArray  *plhs[], 
                 int nrhs, const mxArray  *prhs[] )

{    double   *blksize, *x, *z; 
     int      k, l, blkdim, numblk, cols, idx;

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs != 2){
      mexErrMsgTxt("mexexpand: requires 2 input arguments."); }
   if (nlhs > 1){ 
      mexErrMsgTxt("mexexpand: requires 1 output argument."); }

/* CHECK THE DIMENSIONS */

    numblk = MAX(mxGetM(prhs[0]),mxGetN(prhs[0])); 
    blksize = mxGetPr(prhs[0]);
    cols = 0; 
    for (k=0; k<numblk; k++) { 
        cols = cols + (int)blksize[k]; 
    } 
    x = mxGetPr(prhs[1]); 
    if (mxIsSparse(prhs[1])) { 
        mexErrMsgTxt("mexexpand: sparse x not allowed"); }
    if (MAX(mxGetM(prhs[1]),mxGetN(prhs[1])) != numblk) {
        mexErrMsgTxt("mexexpand: size of blk and x incompatible."); }
    plhs[0] = mxCreateDoubleMatrix(cols,1,mxREAL); 
    z = mxGetPr(plhs[0]);    

    idx = 0; 
    for (k=0; k<numblk; k++) { 
       blkdim = (int)blksize[k]; 
       for (l=0; l<blkdim; l++) { z[idx] = x[k]; idx++; }
    }
    return;
 }
/**********************************************************/

