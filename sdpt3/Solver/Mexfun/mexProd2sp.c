/******************************************************************
* mexProd2sp: compute A*B, where A*B are sparse
* 
*   P = mexProd2sp(A,B,nnz), where P = A*B. 
*   
*   output: P is sparse.
*   
* SDPT3: version 3.0 
* Copyright (c) 1997 by
* K.C. Toh, M.J. Todd, R.H. Tutuncu
* Last Modified: 2 Feb 01
******************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

/**********************************************************
*  A, B sparse.
*  A = mxn, B = nxp  
**********************************************************/
void  prodsp(int m, int n, int p,  
            double *A, mwIndex *irA, mwIndex *jcA,
            double *B, mwIndex *irB, mwIndex *jcB,
            double *P, mwIndex *irP, mwIndex *jcP, double *Btmp)  

{ int  j, k, r, t, Astart, Aend, Bstart, Bend, idx, count;
  double  tmp;  

     count = 0;  jcP[0] = 0; 

     for (j=0; j<p; j++) {
         Bstart = jcB[j];  Bend = jcB[j+1];
         for (k=Bstart; k<Bend; k++) { 
	     tmp = B[k]; 
             r = irB[k]; 
             Astart = jcA[r];  Aend = jcA[r+1]; 
             for (t=Astart; t<Aend; t++) {
                idx = irA[t]; 
                Btmp[idx] += A[t]*tmp; 
	     }
	 }
         for (k=0; k<m; k++) { 
  	     tmp = Btmp[k]; 
             if (tmp != 0) {
	        P[count] = tmp; 
	        irP[count] = k; 
                count++; 
                Btmp[k] = 0; 
	     }
         }
         jcP[j+1]=count;   
     }
     return;
}
/**********************************************************
* 
**********************************************************/
void mexFunction(int nlhs,  mxArray        *plhs[],
           	 int nrhs,  const mxArray  *prhs[] )
{
   double   *A,    *B,    *P; 
   mwIndex  *irA,  *jcA,  *irB,  *jcB, *irP, *jcP;
   double   *Btmp; 
   int       nnz, m1, n1, m2, n2, mlist, nlist, isspA, isspB, k;

/* Check for proper number of arguments */
   if (nrhs != 3){
      mexErrMsgTxt("mexProd2sp: requires 3 input arguments."); }
   else if (nlhs > 2){ 
      mexErrMsgTxt("mexProd2sp: requires 1 output argument."); }
   
/***********************************************/     

       A = mxGetPr(prhs[0]); 
       isspA = mxIsSparse(prhs[0]);
       m1 = mxGetM(prhs[0]); 
       n1 = mxGetN(prhs[0]); 
       if (!isspA) { 
          mexErrMsgTxt("mexProd2sp: A must be sparse"); }        
       if (isspA) { irA = mxGetIr(prhs[0]);
                    jcA = mxGetJc(prhs[0]); }
       B = mxGetPr(prhs[1]); 
       isspB  = mxIsSparse(prhs[1]);
       m2 = mxGetM(prhs[1]); 
       n2 = mxGetN(prhs[1]); 
       if (!isspB) {
          mexErrMsgTxt("mexProd2sp: B must be sparse"); } 
       if (isspB) { irB = mxGetIr(prhs[1]);
                    jcB = mxGetJc(prhs[1]); }
       if (n1!=m2) { 
          mexErrMsgTxt("mexProd2sp: A, B are not compatible"); } 
       nnz = (int) *mxGetPr(prhs[2]); 
       plhs[0] = mxCreateSparse(m1,n2,nnz,mxREAL); 
       P = mxGetPr(plhs[0]); 
       irP = mxGetIr(plhs[0]); 
       jcP = mxGetJc(plhs[0]);
       /**********************************************
        * Do the actual computations in a subroutine 
        **********************************************/

       Btmp = mxCalloc(m1,sizeof(double)); 
       prodsp(m1,n1,n2,A,irA,jcA,B,irB,jcB,P,irP,jcP,Btmp);
       if (isspB) { mxFree(Btmp); }
return;
}		 
/**********************************************************/

