/***********************************************************************
* mexsmat.c : C mex file 
*
*  M = mexsmat(blk,x,isspM,rowidx,colidx,iscellM); 
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
* form Q using the upper triangular part
**********************************************************/
void sym(double *Q, int n)
{
   int j, k, jn, kn; 

   for (j=0; j<n; ++j){
       jn = j*n; 
       for (k=0; k<j ; ++k){ 
           kn = k*n;
           Q[j+kn] = Q[k+jn]; }
   }
   return;
}
/**********************************************************
* Complex case:
* form Q using the upper triangular part
**********************************************************/
void symcmp(double *Q, double *QI, int n)
{
   int j, k, jn, kn; 

   for (j=0; j<n; ++j){
       jn = j*n; 
       for (k=0; k<j ; ++k){ 
           kn = k*n;
           Q[j+kn]  =  Q[k+jn]; 
           QI[j+kn] = -QI[k+jn]; 
       }
   }
   return;
}
/**********************************************************
* form upper triangular part of P
* single block 
**********************************************************/
void smat1(int n, const double ir2, 
           const double *A, const mwIndex *irA, const mwIndex *jcA, int isspA, 
           int mA, int colidx, 
           double *B, mwIndex *irB, mwIndex *jcB, int isspB)

{  int idx, i, j, r, jn, k, kstart, kend, idxj, j2, count;
   double tmp;  
   double hf=0.5; 
   
   if (!isspA & !isspB) { 
      idx = colidx*mA; 
      for (j=0; j<n; j++) { 
          jn = j*n; 
          for (i=0; i<j; i++) { 
              B[i+jn] = ir2*A[idx]; 
              idx++; } 
          B[j+jn] = A[idx];
          idx++; 
      }
   } else if (isspA & !isspB) {      
      j2 = 0; idxj = 0; 
      kstart = jcA[colidx];  kend = jcA[colidx+1]; 
      for (k=kstart; k<kend; k++) { 
          r = irA[k];
          for (j=j2; j<n; j++) {i=r-idxj; if (i>j) {idxj+=j+1;} else {break;}} j2=j; 
          if (i < j) { B[i+j*n] = ir2*A[k]; }
          else       { B[i+j*n] = A[k]; }
      }
   } else if (!isspA & isspB) { 
      idx = colidx*mA; 
      count = 0; 
      for (j=0; j<n; j++) { 
          for (i=0; i<j; i++) {
   	      tmp = A[idx];
              if (tmp != 0) { irB[count] = i; B[count] = ir2*tmp; count++; }
              idx++; 
          }     
          tmp = A[idx];
          if (tmp != 0) { irB[count] = j; B[count] = hf*tmp; count++; }
          idx++; 
          jcB[j+1] = count; 
      }   
   } else if (isspA & isspB) { 
      count = 0; 
      j2 = 0; idxj = 0; 
      kstart = jcA[colidx];  kend = jcA[colidx+1]; 
      for (k=kstart; k<kend; k++) { 
          r = irA[k];
          for (j=j2; j<n; j++) {i=r-idxj; if (i>j) {idxj+=j+1;} else {break;}} j2=j; 
          irB[count] = i;
          if (i<j) {B[count] = ir2*A[k];}  else {B[count] = hf*A[k];}
          ++jcB[j+1]; 
          count++; 
      }   
      for (j=0; j<n; j++) { jcB[j+1] += jcB[j]; }
   }
   if (!isspB) { sym(B,n); }
return; 
}
/**********************************************************
* form upper triangular part of P
* multiple sub-blocks 
**********************************************************/
void smat2(int n, const int numblk, const int *cumblksize, const int *blknnz, 
           const double ir2, 
           const double *A, const mwIndex *irA, const mwIndex *jcA, int isspA, 
           int mA, int colidx, 
           double *B, mwIndex *irB, mwIndex *jcB, int isspB)

{  int idx, i, j, r, jn, k, kstart, kend, idxj, j2, count;
   int t, t2, istart, jstart, jend, rowidx; 
   double tmp;  
   double hf=0.5; 

   if (!isspA) { 
      idx = 0; 
      jstart = 0; jend = 0; 
      for (t=0; t<numblk; t++) {   	
 	  jend = cumblksize[t+1]; 
          istart = jstart;
          idxj = colidx*mA; 
          for (j=jstart; j<jend; j++){
   	      idxj += j-jstart; 
              rowidx = blknnz[t]-istart+idxj; 
              for (i=istart; i<j; i++) {
  		  irB[idx] = i;  B[idx] = ir2*A[rowidx+i]; idx++; }                    
              irB[idx] = j;  B[idx] = hf*A[rowidx+j]; idx++; 
              jcB[j+1] = idx; 
          } 
          jstart = jend; 
      }
   } else {
      jstart = 0; jend = cumblksize[1]; t2 = 0; 
      kstart = jcA[colidx]; kend = jcA[colidx+1]; 
      count  = 0; j2 = 0; idxj = 0; 
      for (k=kstart; k<kend; k++) { 
          r = irA[k];
          for (t=t2; t<numblk; t++) { if (r-blknnz[t+1]<0) {break;} } 
          if (t > t2) { 
             t2 = t; 
             jstart = cumblksize[t2]; jend = cumblksize[t2+1]; 
             idxj = blknnz[t2];
             j2 = jstart; 
          }
          for (j=j2; j<jend; j++) {i=jstart+r-idxj; if (i>j) {idxj+=j-jstart+1;} else {break;}}
          j2=j; 
          irB[count] = i;
          if (i<j) {B[count] = ir2*A[k];}  else {B[count] = hf*A[k];}
          ++jcB[j+1]; 
          count++; 
      }  
      for (j=0; j<n; j++) { jcB[j+1] += jcB[j]; }
   }
return; 
}
/**********************************************************
* Complex case:
* form upper triangular part of P
* single block 
**********************************************************/
void smat1cmp(int n, const double ir2, 
           const double *A, const mwIndex *irA, const mwIndex *jcA, int isspA, 
           int mA, int colidx, 
           double *B, mwIndex *irB, mwIndex *jcB, int isspB,
           double *AI, double *BI)

{  int idx, i, j, r, jn, k, kstart, kend, idxj, j2, count, ind;
   double tmp, tmp2;  
   double hf=0.5; 
   
   if (!isspA & !isspB) { 
      idx = colidx*mA; 
      for (j=0; j<n; j++) { 
          jn = j*n; 
          for (i=0; i<j; i++) { 
              B[i+jn]  = ir2*A[idx]; 
              BI[i+jn] = ir2*AI[idx]; 
              idx++; } 
          B[j+jn]  = A[idx];
          BI[j+jn] = AI[idx];
          idx++; 
      }
   } else if (isspA & !isspB) {      
      j2 = 0; idxj = 0; 
      kstart = jcA[colidx];  kend = jcA[colidx+1]; 
      for (k=kstart; k<kend; k++) { 
          r = irA[k];
          for (j=j2; j<n; j++) {i=r-idxj; if (i>j) {idxj+=j+1;} else {break;}} j2=j; 
          if (i < j) { ind = i+j*n; B[ind] = ir2*A[k]; BI[ind] = ir2*AI[k]; }
          else       { ind = i+j*n; B[ind] = A[k];     BI[ind] = AI[k];}
      }
   } else if (!isspA & isspB) { 
      idx = colidx*mA; 
      count = 0; 
      for (j=0; j<n; j++) { 
          for (i=0; i<j; i++) {
   	      tmp = A[idx]; tmp2 = AI[idx]; 
              if ((tmp != 0) || (tmp2 != 0)) { 
                 irB[count] = i; 
                 B[count] = ir2*tmp; BI[count] = ir2*tmp2; 
                 count++; 
	      }
              idx++; 
          }     
          tmp = A[idx]; tmp2 = AI[idx]; 
          if ((tmp != 0) || (tmp2 != 0)) { 
             irB[count] = j; B[count] = hf*tmp; BI[count] = hf*tmp2; 
             count++; 
          }
          idx++; 
          jcB[j+1] = count; 
      }   
   } else if (isspA & isspB) { 
      count = 0; 
      j2 = 0; idxj = 0; 
      kstart = jcA[colidx];  kend = jcA[colidx+1]; 
      for (k=kstart; k<kend; k++) { 
          r = irA[k];
          for (j=j2; j<n; j++) {i=r-idxj; if (i>j) {idxj+=j+1;} else {break;}} j2=j; 
          irB[count] = i;
          if (i<j) {
             B[count] = ir2*A[k]; BI[count] = ir2*AI[k];
          }  else {
             B[count] = hf*A[k];  BI[count] = hf*AI[k];
	  }
          ++jcB[j+1]; 
          count++; 
      }   
      for (j=0; j<n; j++) { jcB[j+1] += jcB[j]; }
   }
   if (!isspB) { symcmp(B,BI,n); }
return; 
}
/**********************************************************
* Complex case:
* form upper triangular part of P
* multiple sub-blocks 
**********************************************************/
void smat2cmp(int n, const int numblk, const int *cumblksize, const int *blknnz, 
           const double ir2, 
           const double *A, const mwIndex *irA, const mwIndex *jcA, int isspA, 
           int mA, int colidx, 
           double *B, mwIndex *irB, mwIndex *jcB, int isspB,
           double *AI, double *BI)

{  int  idx, i, j, r, jn, k, kstart, kend, idxj, j2, count;
   int  t, t2, istart, jstart, jend, rowidx; 
   double tmp;  
   double hf=0.5; 

   if (!isspA) { 
      idx = 0; 
      jstart = 0; jend = 0; 
      for (t=0; t<numblk; t++) {   	
 	  jend = cumblksize[t+1]; 
          istart = jstart;
          idxj = colidx*mA; 
          for (j=jstart; j<jend; j++){
   	      idxj += j-jstart; 
              rowidx = blknnz[t]-istart+idxj; 
              for (i=istart; i<j; i++) {
  		  irB[idx] = i;  
                  B[idx] = ir2*A[rowidx+i];  BI[idx] = ir2*AI[rowidx+i]; 
                  idx++; 
              }                    
              irB[idx] = j;  
              B[idx] = hf*A[rowidx+j];  BI[idx] = hf*AI[rowidx+j]; 
              idx++; 
              jcB[j+1] = idx; 
          } 
          jstart = jend; 
      }
   } else {
      jstart = 0; jend = cumblksize[1]; t2 = 0; 
      kstart = jcA[colidx]; kend = jcA[colidx+1]; 
      count  = 0; j2 = 0; idxj = 0; 
      for (k=kstart; k<kend; k++) { 
          r = irA[k];
          for (t=t2; t<numblk; t++) { if (r-blknnz[t+1]<0) {break;} } 
          if (t > t2) { 
             t2 = t; 
             jstart = cumblksize[t2]; jend = cumblksize[t2+1]; 
             idxj = blknnz[t2];
             j2 = jstart; 
          }
          for (j=j2; j<jend; j++) {i=jstart+r-idxj; if (i>j) {idxj+=j-jstart+1;} else {break;}}
          j2=j; 
          irB[count] = i;
          if (i<j) {
             B[count] = ir2*A[k]; BI[count] = ir2*AI[k];
          }  else {
             B[count] = hf*A[k];  BI[count] = hf*AI[k];
          }
          ++jcB[j+1]; 
          count++; 
      }  
      for (j=0; j<n; j++) { jcB[j+1] += jcB[j]; }
   }
return; 
}
/**********************************************************
* 
***********************************************************/
void mexFunction(int nlhs,   mxArray  *plhs[], 
                 int nrhs,   const mxArray  *prhs[] )

{    mxArray  *rhs[2]; 
     mxArray  *blk_cell_pr, *A_cell_pr;
     double   *A,  *B,  *blksize, *AI, *BI;
     mwIndex  *irA, *jcA, *irB, *jcB;
     int      *cumblksize, *blknnz;
     int       mblk, mA, nA, m1, n1, rowidx, colidx, isspA, isspB;
     int       iscellA, iscmpA; 

     mwIndex  subs[2];
     mwSize   nsubs=2; 
     int      n, n2, k, nsub, index, numblk, NZmax;
     double   ir2=1/sqrt(2); 

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

     if (nrhs < 2){
         mexErrMsgTxt("mexsmat: requires at least 2 input arguments."); }
     else if (nlhs>1){ 
         mexErrMsgTxt("mexsmat: requires 1 output argument."); }

/* CHECK THE DIMENSIONS */

     iscellA = mxIsCell(prhs[1]); 
     if (iscellA) { mA = mxGetM(prhs[1]); nA = mxGetN(prhs[1]); }
     else         { mA = 1; nA = 1; }
     if (mxGetM(prhs[0]) != mA) {
         mexErrMsgTxt("mexsmat: blk and Avec not compatible"); }

/***** main body *****/ 
       
     if (nrhs > 3) {rowidx = (int)*mxGetPr(prhs[3]); } else {rowidx = 1;}  
     if (rowidx > mA) {
         mexErrMsgTxt("mexsmat: rowidx exceeds size(Avec,1)."); }
     subs[0] = rowidx-1;  /* subtract 1 to adjust for Matlab index */
     subs[1] = 1;
     index = mxCalcSingleSubscript(prhs[0],nsubs,subs); 
     blk_cell_pr = mxGetCell(prhs[0],index);
     if (blk_cell_pr == NULL) { 
        mexErrMsgTxt("mexsmat: blk not properly specified"); }    
     numblk  = mxGetN(blk_cell_pr);            
     blksize = mxGetPr(blk_cell_pr);
     cumblksize = mxCalloc(numblk+1,sizeof(int)); 
     blknnz = mxCalloc(numblk+1,sizeof(int)); 
     cumblksize[0] = 0; blknnz[0] = 0; 
     n = 0;  n2 = 0; 
     for (k=0; k<numblk; ++k) {
          nsub = (int) blksize[k];
          n  += nsub; 
          n2 += nsub*(nsub+1)/2;  
          cumblksize[k+1] = n; 
          blknnz[k+1] = n2;
     }
     /***** assign pointers *****/
     if (iscellA) { 
         subs[0] = rowidx-1; 
         subs[1] = 0;
         index = mxCalcSingleSubscript(prhs[1],nsubs,subs); 
         A_cell_pr = mxGetCell(prhs[1],index); 
         A  = mxGetPr(A_cell_pr); 
         m1 = mxGetM(A_cell_pr); 
         n1 = mxGetN(A_cell_pr);
         isspA = mxIsSparse(A_cell_pr);
         iscmpA = mxIsComplex(A_cell_pr);  
         if (isspA) { irA = mxGetIr(A_cell_pr);
                      jcA = mxGetJc(A_cell_pr); } 
         if (iscmpA) { AI = mxGetPi(A_cell_pr); }
     } else { 
         A  = mxGetPr(prhs[1]); 
         m1 = mxGetM(prhs[1]); 
         n1 = mxGetN(prhs[1]); 
         isspA = mxIsSparse(prhs[1]); 
         iscmpA = mxIsComplex(prhs[1]);  
         if (isspA) {  irA = mxGetIr(prhs[1]); 
                       jcA = mxGetJc(prhs[1]); }
         if (iscmpA) { AI = mxGetPi(prhs[1]); }
     }
     if (numblk > 1) { 
        isspB = 1; 
     } else { 
        if (nrhs > 2) {isspB = (int)*mxGetPr(prhs[2]);} else {isspB = isspA;} 
     }
     if (nrhs > 4) {colidx = (int)*mxGetPr(prhs[4]) -1;} else {colidx = 0;} 
     if (colidx > n1) { 
         mexErrMsgTxt("mexsmat: colidx exceeds size(Avec,2)."); 
     }    
     /***** create return argument *****/
     if (isspB) {
	if (isspA) { NZmax = jcA[colidx+1]-jcA[colidx]; } 
        else       { NZmax = blknnz[numblk]; } 
        if (iscmpA) { 
	   rhs[0] = mxCreateSparse(n,n,2*NZmax,mxCOMPLEX); 
        } else {
           rhs[0] = mxCreateSparse(n,n,2*NZmax,mxREAL); 
        }
	B = mxGetPr(rhs[0]); 
        irB = mxGetIr(rhs[0]); 
        jcB = mxGetJc(rhs[0]);
        if (iscmpA) { BI = mxGetPi(rhs[0]); }
     } else {
        if (iscmpA) { 
          plhs[0] = mxCreateDoubleMatrix(n,n,mxCOMPLEX); 
        } else {
          plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL); 
        }
        B = mxGetPr(plhs[0]); 
        if (iscmpA) { BI = mxGetPi(plhs[0]); }
     }
     /***** Do the computations in a subroutine *****/
     if (iscmpA) { 
       if (numblk == 1) { 
          smat1cmp(n,ir2,A,irA,jcA,isspA,m1,colidx,B,irB,jcB,isspB,AI,BI);  
       } else { 
          smat2cmp(n,numblk,cumblksize,blknnz,ir2,A,irA,jcA,isspA,m1,colidx,B,irB,jcB,isspB,AI,BI);
       }
     } else {
       if (numblk == 1) { 
          smat1(n,ir2,A,irA,jcA,isspA,m1,colidx,B,irB,jcB,isspB);  
       } else { 
          smat2(n,numblk,cumblksize,blknnz,ir2,A,irA,jcA,isspA,m1,colidx,B,irB,jcB,isspB);
       }
     }
     if (isspB) {
        /*** if isspB, (actual B) = B+B' ****/ 
        mexCallMATLAB(1, &rhs[1], 1, &rhs[0], "ctranspose"); 
        mexCallMATLAB(1, &plhs[0],2, rhs, "+");  
        mxDestroyArray(*rhs); 
     }
     mxFree(blknnz); 
     mxFree(cumblksize);
 return;
 }
/**********************************************************/

