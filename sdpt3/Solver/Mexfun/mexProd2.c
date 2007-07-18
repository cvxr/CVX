/******************************************************************
* mexProd2.c : C-MEX file to compute the product of two matrices. 
*       
*         P = mexProd2(blk,A,B,type)
*
* input:  blk = 1x2 cell array describing the block structure of A and B
*         A = mxn matrix.
*         B = nxp matrix.
*      type = 0     general matrix product
*             1     if P is symmetric
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
* saxpy:  z = z + alpha*y
**********************************************************/
void saxpy(double x, double *y, int idx1, 
           double *z, int idx2, int istart, int iend)
{  int i;
  
   for(i=istart; i<iend-3; i++){             /* LEVEL 4 */
      z[i+idx2] += x * y[i+idx1]; i++;
      z[i+idx2] += x * y[i+idx1]; i++;
      z[i+idx2] += x * y[i+idx1]; i++;
      z[i+idx2] += x * y[i+idx1]; 
   }
   if(i < iend-1){                           /* LEVEL 2 */
      z[i+idx2] += x * y[i+idx1]; i++;
      z[i+idx2] += x * y[i+idx1]; i++;
   }
   if(i < iend){                             /* LEVEL 1 */
      z[i+idx2] += x * y[i+idx1]; 
   }
   return;
}
/**********************************************************
* form P using the upper triangular part
**********************************************************/
void symmetrize(double *P, int n)

{  int j, k, jn; 

   for (j=0; j<n; j++){
       jn = j*n; 
       for (k=0; k<j; k++){ P[j+k*n] = P[k+jn]; }
   }
   return;
}
/**********************************************************
* A dense, B dense
**********************************************************/
void product(double *A, double *B, double *P, 
             int m, int n, int p, int type)

{  int  i, j, k, jm, jn, km, kstart, kend;
   int     istart, iend;
   double  tmp; 

      for (j=0; j<p; j++){
          kstart = 0;  kend = n;
          jm = j*m;  jn = j*n; 
	  for (k=kstart; k<kend; k++){
              istart = 0; 
              if (type==1) {iend = j+1;} else {iend = m;} 
              tmp = B[k+jn];
              if (tmp != 0) {
                 km = k*m;  
                 saxpy(tmp,A,km,P,jm,istart,iend); }	    
	  }
      }
      if (type==1) { symmetrize(P,m); }
      return; 
}
/**********************************************************
* A dense, B sparse
**********************************************************/
void product2(double *A, double *B, mwIndex *irB, mwIndex *jcB, 
              double *P, int m, int n, int p, int type)

{  int  i, j, k, r, kstart, kend, istart, iend, jm, rm;
   double  tmp;
    
      for (j=0; j<p; j++){
          kstart = jcB[j]; 
          kend   = jcB[j+1]; 
          jm = j*m; 
	  for (k=kstart; k<kend; k++){
              r = irB[k]; 
              tmp = B[k];
              istart = 0;
              if (type==1) {iend = j+1;} else {iend = m;}
              if (tmp != 0) {
                 rm = r*m; 
		 saxpy(tmp,A,rm,P,jm,istart,iend); }
          }
      }
      if (type==1) { symmetrize(P,m); }
      return;
}
/**********************************************************
* A sparse, B dense
**********************************************************/
void product3(double *A, mwIndex *irA, mwIndex *jcA, double *B, 
              double *P, int m, int n, int p, int type)
                     
{  int  i, j, k, ri, kstart, kend, istart, iend, jm, jn, sym;
   double  tmp;
    
      for (j=0; j<p; j++){
          kstart = 0;  kend = n;  
          if (type==1) {sym = 1;} else {sym = 0;} 
          jm = j*m;  jn = j*n; 
	  for (k=kstart; k<kend; k++){
              tmp  = B[k+jn];
              istart = jcA[k]; 
              iend = jcA[k+1]; 
              if (tmp != 0) {
  		 for (i=istart; i<iend; i++) {
                    ri = irA[i]; 
                    if (ri > j & sym) { break; }
		    P[ri+jm] += tmp*A[i]; }
	      }
	  }
      }
      if (type==1) { symmetrize(P,m); }
      return;
}
/**********************************************************
* A sparse, B sparse
**********************************************************/
void product4(double *A, mwIndex *irA, mwIndex *jcA, 
              double *B, mwIndex *irB, mwIndex *jcB, 
              double *P, mwIndex *irP, mwIndex *jcP,
              double *Ptmp, int numblk, int *cumblk)

{  int  i, j, k, l, r, t, istart, iend, kstart, kend, jstart, jend;
   int  idx;
   double  tmp; 

     idx = 0;  jcP[0]=0;
     for (l=0; l<numblk; l++) { 
        jstart = cumblk[l]; jend = cumblk[l+1];
        for (j=jstart; j<jend; j++){
            kstart = jcB[j];  kend  = jcB[j+1];
           /**** forming jth column of P ****/
           for (k=kstart; k<kend; k++) {
               r = irB[k];
               tmp = B[k];
               istart = jcA[r]; iend = jcA[r+1];
               for (i=istart; i<iend; i++) {
                    t = irA[i];
	            Ptmp[t] += tmp*A[i]; }
	   }
	   for (k=jstart; k<jend; k++) {
                tmp = Ptmp[k]; 
	        if (tmp != 0) {
                    P[idx] = tmp; irP[idx] = k;
                    Ptmp[k] = 0;  idx++; }
	   } 
           jcP[j+1] = idx;    
	}
     }
     jcP[jend] = idx; 
     return; 
}
/**********************************************************
* elementwise product of two real column vectors. 
**********************************************************/
void product5(double *A, mwIndex *irA, mwIndex *jcA, 
              double *B, mwIndex *irB, mwIndex *jcB, double *P, 
              int n, int isspA, int isspB)

{  int  k, kx, ky, kx2, ky2, rx, ry;

      if ( !isspA & !isspB ) {
         for (k=0; k<n; k++){ P[k] = A[k]*B[k]; }
      }
      else if ( isspA & !isspB ) {
         kx = jcA[0];  kx2 = jcA[1];
	 for (k=kx; k<kx2; k++) {
             rx = irA[k]; 
             P[rx] = A[k]*B[rx]; }   
      }
      else if ( !isspA & isspB ) {
         ky = jcB[0];  ky2 = jcB[1];
	 for (k=ky; k<ky2; k++) {
             ry = irB[k]; 
             P[ry] = A[ry]*B[k]; }   
      }
      else if ( isspA & isspB ) {
           kx = jcA[0];  kx2 = jcA[1];  rx = irA[kx];
           ky = jcB[0];  ky2 = jcB[1];  ry = irB[ky]; 
           while ( (kx<kx2) & (ky<ky2) ){
               if (rx == ry) { 
		  P[rx] = A[kx]*B[ky]; 
                  kx++; ky++; 
                  rx = irA[kx];
                  ry = irB[ky]; } 
               else if (rx < ry) { 
                  kx++; 
                  rx = irA[kx]; }
               else { 
                  ky++; 
                  ry = irB[ky]; }                   
	   }
      }
      return; 
}
/**********************************************************/
void mexFunction(int nlhs,  mxArray        *plhs[],
	         int nrhs,  const mxArray  *prhs[] )
{
   mxArray  *blk_cell_pr;
   double   *A,   *B,   *P, *blksize, *Ptmp; 
   mwIndex  *irA, *jcA, *irB, *jcB, *irP, *jcP;
   int      *cumblk;
   int       isspA, isspB, m1, n1, m2, n2;
   int       type, index, numblk, NZmax, cols, i, l;
   mwIndex   subs[2];
   mwSize    nsubs=2;

/* Check for proper number of arguments */
   if (nrhs<3){
      mexErrMsgTxt("mexProd2: requires at least 3 input arguments."); }
   else if (nlhs>2){ 
      mexErrMsgTxt("mexProd2: requires 1 output argument."); }
   if (mxIsCell(prhs[1]) || mxIsCell(prhs[2])) {
      mexErrMsgTxt("mexProd2: 2ND and 3RD input must both be matrices"); }
   if (mxGetM(prhs[0]) > 1) { 
      mexErrMsgTxt("mexProd2: blk can only have 1 row"); }     

/*** get pointers ***/  

    if (nrhs > 3) { type = (int)mxGetScalar(prhs[3]); }
    else          { type = 0; }
    subs[0] = 0;  subs[1] = 1;
    index = mxCalcSingleSubscript(prhs[0],nsubs,subs); 
    blk_cell_pr = mxGetCell(prhs[0],index); 
    blksize = mxGetPr(blk_cell_pr); 
    numblk = mxGetN(blk_cell_pr);
    cumblk = mxCalloc(numblk+1,sizeof(int)); 
    NZmax = 0; 
    for (l=0; l<numblk; l++) { 
        cols = (int)blksize[l]; 
        cumblk[l+1] = cumblk[l] + cols;  
        NZmax += cols*cols; }
    A = mxGetPr(prhs[1]); 
    m1 = mxGetM(prhs[1]); 
    n1 = mxGetN(prhs[1]); 
    isspA = mxIsSparse(prhs[1]);
    if (isspA) { irA = mxGetIr(prhs[1]);
                 jcA = mxGetJc(prhs[1]); }
    B = mxGetPr(prhs[2]); 
    m2 = mxGetM(prhs[2]); 
    n2 = mxGetN(prhs[2]); 
    isspB = mxIsSparse(prhs[2]);
    if (isspB) { irB = mxGetIr(prhs[2]);
        	 jcB = mxGetJc(prhs[2]); }       
    if ((n1!=m2) & !(n1==1 & n2==1)) { 
        mexErrMsgTxt("mexProd2: 2ND and 3RD input not compatible"); }
    if ((numblk > 1) & !(isspA & isspB) & !(n1==1 & n2==1)) { 
       mexErrMsgTxt("mexProd2: 2ND and 3RD must be both sparse"); }

/***** create return argument *****/   

    if (isspA & isspB & !(n1==1 & n2==1)){ 
       plhs[0] = mxCreateSparse(m1,n2,NZmax,mxREAL); 
       P = mxGetPr(plhs[0]); irP = mxGetIr(plhs[0]); jcP = mxGetJc(plhs[0]); }
    else { 
       plhs[0] = mxCreateDoubleMatrix(m1,n2,mxREAL);
       P = mxGetPr(plhs[0]); 
    }
    if (isspA & isspB & !(n1==1 & n2==1)) {  
       Ptmp = mxCalloc(cumblk[numblk],sizeof(double)); 
    }
/**********************************************
* Do the actual computations in a subroutine 
**********************************************/

    if (m1 == m2 & n1 == 1 & n2 == 1) {
       product5(A, irA, jcA, B, irB, jcB, P, m1, isspA, isspB); 
    } else {
      if (!isspA & !isspB){
         product(A, B, P, m1, n1, n2, type); }
      else if (!isspA & isspB){ 
         product2(A, B, irB, jcB, P, m1, n1, n2, type); }
      else if (isspA & !isspB){ 
         product3(A, irA, jcA, B, P, m1, n1, n2, type); }
      else if (isspA & isspB){
         product4(A, irA, jcA, B, irB, jcB,P,irP,jcP,Ptmp,numblk,cumblk); 
      } 
    }
   mxFree(cumblk); 
   if (isspA & isspB & !(n1==1 & n2==1)) { mxFree(Ptmp); }
   return;
 }
/**********************************************************/

