/***********************************************************************
* mexinprod.c : C mex file to compute inner product of 2 vectors,  
*           
* Ax = mexinprod(blk,Avec,x,k,p); 
* where Ax = Avec{p}(:,1:k)'*x 
*
* SDPT3: version 3.0
* Copyright (c) 1997 by
* K.C. Toh, M.J. Todd, R.H. Tutuncu
* Last Modified: 2 Feb 01   
************************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

/**********************************************************
* realdot1: x dense, y dense 
**********************************************************/
double realdot1(const double *x, const int col, 
                const double *y, const int n)

{  int i, idx;
   double r; 

   idx=col*n; 
   r=0.0;
   for (i=0; i<n-3; i++) {             /* LEVEL 4 */
      r += x[i+idx] * y[i]; i++;
      r += x[i+idx] * y[i]; i++;
      r += x[i+idx] * y[i]; i++;
      r += x[i+idx] * y[i]; }
   if (i<n-1) {                        /* LEVEL 2 */
      r += x[i+idx] * y[i]; i++;
      r += x[i+idx] * y[i]; i++; }
   if (i<n) {                          /* LEVEL 1 */
      r += x[i+idx] * y[i]; }
   return r; 
}
/**********************************************************
* realdot2: x sparse, y dense 
**********************************************************/
double realdot2(const double *x, mwIndex *irx, mwIndex *jcx, 
                const int col, const double *y)

{  int i, row, idxstart, idxend;
   double r; 

   idxstart=jcx[col]; idxend=jcx[col+1]; 
   r=0.0;
   for (i=idxstart; i<idxend-3; i++) {      /* LEVEL 4 */
      row = irx[i]; r += x[i]*y[row]; i++;
      row = irx[i]; r += x[i]*y[row]; i++;
      row = irx[i]; r += x[i]*y[row]; i++;
      row = irx[i]; r += x[i]*y[row]; }
   if (i<idxend-1) {                        /* LEVEL 2 */
      row = irx[i]; r += x[i]*y[row]; i++;
      row = irx[i]; r += x[i]*y[row]; i++; }
   if (i<idxend) {                          /* LEVEL 1 */
      row = irx[i]; r += x[i]*y[row]; }
   return r; 
}
/**********************************************************/
void mexFunction(int nlhs,   mxArray  *plhs[], 
                 int nrhs,   const mxArray  *prhs[] )

{    mxArray  *A_cell_pr;
     double   *A, *B;
     mwIndex  *irA, *jcA, *irB, *jcB; 
     double   *Btmp, *tr;
     int       isspA, isspB, iscellA, iscellB;

     mwIndex   subs[2];
     mwSize    nsubs=2; 
     int       mA, nA, m1, n1, m2, n2, j, index;    
     int       rowidx, colidx, r, k, kstart, kend;


/* CHECK THE DIMENSIONS */

   iscellA = mxIsCell(prhs[1]); 
   mA = mxGetM(prhs[1]); 
   if (!iscellA) { mA = 1; }  
   if (nrhs < 2) {
       mexErrMsgTxt(" mexinprod: must have at least 3 inputs"); }
   if (nlhs>2) { 
       mexErrMsgTxt("mexinprod: requires 1 output argument"); }
   if (nrhs > 4) { rowidx = (int)*mxGetPr(prhs[4]); } else { rowidx = 1; }
   if (rowidx > mA) {
      mexErrMsgTxt("mexinprod: rowidx exceeds size(Avec,1)"); 
   }
   /***** assign pointers *****/
    
       iscellB = mxIsCell(prhs[2]); 
       isspB = mxIsSparse(prhs[2]);
       m2 = mxGetM(prhs[2]); 
       n2 = mxGetN(prhs[2]);            
       if ((n2 > 1) || (iscellB)) { 
 	  mexErrMsgTxt("mexinprod: 3RD input must be a column vector"); }
       if (isspB) { 
          irB = mxGetIr(prhs[2]);
          jcB = mxGetJc(prhs[2]); 
          Btmp = mxGetPr(prhs[2]); 
          /***** copy Btmp to B *****/ 
          B = mxCalloc(m2,sizeof(double)); 
          kstart = jcB[0]; kend = jcB[1]; 
          for (k=kstart; k<kend; k++) { 
	      r = irB[k]; B[r] = Btmp[k]; }
       } else {
          B = mxGetPr(prhs[2]); 
       }       
       if (iscellA) {
  	  subs[0] = rowidx-1;  /* subtract 1 to adjust for Matlab index */
          subs[1] = 0;
          index = mxCalcSingleSubscript(prhs[1],nsubs,subs); 
          A_cell_pr = mxGetCell(prhs[1],index); 
          A  = mxGetPr(A_cell_pr); 
          m1 = mxGetM(A_cell_pr); 
          n1 = mxGetN(A_cell_pr);
          isspA = mxIsSparse(A_cell_pr);               
          if (isspA) { irA = mxGetIr(A_cell_pr);
                       jcA = mxGetJc(A_cell_pr); } 
       } else {
          A = mxGetPr(prhs[1]); 
          m1 = mxGetM(prhs[1]); 
          n1 = mxGetN(prhs[1]);
          isspA = mxIsSparse(prhs[1]);
          if (isspA) { irA = mxGetIr(prhs[1]);
	               jcA = mxGetJc(prhs[1]); }
       }           
       if (nrhs > 3) { colidx = (int)*mxGetPr(prhs[3]); } else { colidx = 1; }
       if (colidx > n1) { 
          mexErrMsgTxt("mexinprod: colidx exceeds size(Avec,2)"); }
       if (m1 != m2) {
          mexErrMsgTxt("mexinprod: 2ND and 3RD input not compatible."); 
       }
       /***** create return argument *****/

       plhs[0] = mxCreateDoubleMatrix(colidx,1,mxREAL); 
       tr = mxGetPr(plhs[0]); 

       /***** compute <Aj,B> *****/
       if (isspA) { 
          for (j=0; j<colidx; j++){
	      tr[j] = realdot2(A,irA,jcA,j,B); }      
       } else {
          for (j=0; j<colidx; j++){
	      tr[j] = realdot1(A,j,B,m1); }       
       }
       if (isspB) { mxFree(B); } 
 return;
}
/**********************************************************/

