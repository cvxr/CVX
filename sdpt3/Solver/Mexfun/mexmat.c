/***********************************************************************
* mexmat.c : C mex file 
*
*  M = mexmat(blk,x,isspM,rowidx,colidx,iscellM); 
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
* single block (real)
**********************************************************/
void  mat1(int n, 
           double *A, mwIndex *irA, mwIndex *jcA, int isspA,
           int mA, int colidx, 
           double *B, mwIndex *irB, mwIndex *jcB, int isspB)

{  int idx, i, j, r, jn, k, kstart, kend, idxj, j2, count;
   double tmp;  
   
   if (!isspA & !isspB) { 
      idx = colidx*mA; 
      for (j=0; j<n; j++) { 
          jn = j*n; 
          for (i=0; i<n; i++) { 
              B[i+jn] = A[idx]; 
              idx++; } 
      }
   } else if (!isspA & isspB) { 
      idx = colidx*mA; 
      count = 0; 
      for (j=0; j<n; j++) { 
          for (i=0; i<n; i++) {
   	      tmp = A[idx];
              if (tmp != 0) { irB[count] = i; B[count] = tmp; count++; }
              idx++; 
          }     
          jcB[j+1] = count; 
      }   
   } else if (isspA & !isspB) {
      j2 = 0; idxj = 0; 
      kstart = jcA[colidx];  kend = jcA[colidx+1]; 
      for (k=kstart; k<kend; k++) { 
          r = irA[k];
          for (j=j2; j<n; j++) {i=r-idxj; if (i>=n) {idxj+=n;} else {break;}} j2=j; 
          B[i+j*n] = A[k]; 
      }
   } else if (isspA & isspB) {
      count = 0; 
      j2 = 0; idxj = 0; 
      kstart = jcA[colidx];  kend = jcA[colidx+1]; 
      for (k=kstart; k<kend; k++) { 
          r = irA[k];
          for (j=j2; j<n; j++) {i=r-idxj; if (i>=n) {idxj+=n;} else {break;}} j2=j; 
          irB[count] = i;
          B[count] = A[k];
          ++jcB[j+1]; 
          count++; 
      }   
      for (j=0; j<n; j++) { jcB[j+1] += jcB[j]; }
   }
return; 
}
/**********************************************************
* single block (complex)
**********************************************************/
void  mat1cmp(int n, 
           double *A, mwIndex *irA, mwIndex *jcA, int isspA,
           int mA, int colidx, 
           double *B, mwIndex *irB, mwIndex *jcB, int isspB,
           double *AI, double *BI)

{  int idx, ind, i, j, r, jn, k, kstart, kend, idxj, j2, count;
   double tmp, tmp2;  
   
   if (!isspA & !isspB) { 
      idx = colidx*mA; 
      for (j=0; j<n; j++) { 
          jn = j*n; 
          for (i=0; i<n; i++) { 
              B[i+jn]  = A[idx]; 
              BI[i+jn] = AI[idx]; 
              idx++; } 
      }
   } else if (!isspA & isspB) { 
      idx = colidx*mA; 
      count = 0; 
      for (j=0; j<n; j++) { 
          for (i=0; i<n; i++) {
	      tmp = A[idx]; tmp2 = AI[idx];
              if ((tmp != 0) || (tmp2 != 0)) { 
                 irB[count] = i; B[count]=tmp; BI[count]=tmp2; count++; 
              }
              idx++; 
          }     
          jcB[j+1] = count; 
      }   
   } else if (isspA & !isspB) {
      j2 = 0; idxj = 0; 
      kstart = jcA[colidx];  kend = jcA[colidx+1]; 
      for (k=kstart; k<kend; k++) { 
          r = irA[k];
          for (j=j2; j<n; j++) {i=r-idxj; if (i>=n) {idxj+=n;} else {break;}} j2=j; 
          ind = i+j*n; 
          B[ind] = A[k]; BI[ind] = AI[k]; 
      }
   } else if (isspA & isspB) {
      count = 0; 
      j2 = 0; idxj = 0; 
      kstart = jcA[colidx];  kend = jcA[colidx+1]; 
      for (k=kstart; k<kend; k++) { 
          r = irA[k];
          for (j=j2; j<n; j++) {i=r-idxj; if (i>=n) {idxj+=n;} else {break;}} j2=j; 
          irB[count] = i;
          B[count] = A[k]; BI[count] = AI[k]; 
          ++jcB[j+1]; 
          count++; 
      }   
      for (j=0; j<n; j++) { jcB[j+1] += jcB[j]; }
   }
return; 
}
/**********************************************************
* B is sparse 
* multiple sub-blocks (real)
**********************************************************/
void  mat2(int n, int numblk, int *cumblksize, int *blknnz, 
           double *A, int mA, int colidx,  
           double *B, mwIndex *irB, mwIndex *jcB, int isspB)

{  int idx, i, j, r, jn, k, kstart, kend, idxj, j2, count;
   int t, t2, istart, jstart, jend, rowidx, nsub; 
   double tmp;  

      idx = 0; 
      jstart = 0; jend = 0; jcB[0]=0;  
      for (t=0; t<numblk; t++) {   	
 	  jend = cumblksize[t+1]; 
          istart = jstart;
          idxj = colidx*mA; 
          nsub = jend-jstart; 
          for (j=jstart; j<jend; j++){
   	      idxj = (j-jstart)*nsub; 
              rowidx = blknnz[t]-istart+idxj; 
              for (i=jstart; i<jend; i++) {
  		  irB[idx] = i;  B[idx] = A[rowidx+i]; idx++; }                    
              jcB[j+1] = idx; 
          } 
          jstart = jend; 
      }
return; 
}
/**********************************************************
* B is sparse 
* multiple sub-blocks (complex)
**********************************************************/
void  mat2cmp(int n, int numblk, int *cumblksize, int *blknnz, 
           double *A, int mA, int colidx,  
           double *B, mwIndex *irB, mwIndex *jcB, int isspB,
           double *AI, double *BI)

{  int idx, i, j, r, jn, k, kstart, kend, idxj, j2, count;
   int t, t2, istart, jstart, jend, rowidx, nsub; 
   double tmp;  

      idx = 0; 
      jstart = 0; jend = 0; jcB[0]=0;  
      for (t=0; t<numblk; t++) {   	
 	  jend = cumblksize[t+1]; 
          istart = jstart;
          idxj = colidx*mA; 
          nsub = jend-jstart; 
          for (j=jstart; j<jend; j++){
   	      idxj = (j-jstart)*nsub; 
              rowidx = blknnz[t]-istart+idxj; 
              for (i=jstart; i<jend; i++) {
  		  irB[idx] = i;  
                  B[idx]  = A[rowidx+i]; 
                  BI[idx] = AI[rowidx+i]; 
                  idx++; 
              }                    
              jcB[j+1] = idx; 
          } 
          jstart = jend; 
      }
return; 
}
/**********************************************************
* 
***********************************************************/
void mexFunction(int nlhs, mxArray  *plhs[], 
                 int nrhs, const mxArray  *prhs[] )

{    
     mxArray  *blk_cell_pr, *A_cell_pr;
     double   *A, *B, *AI, *BI, *blksize, *Atmp, *AItmp;
     mwIndex  *irA, *jcA, *irB, *jcB;
     int      *cumblksize, *blknnz;
     int       iscellA, mblk, mA, nA, m1, n1, rowidx, colidx, isspA, isspB;
     int       iscmpA, iscmpB; 

     mwIndex  subs[2];
     mwSize   nsubs=2; 
     int      n, n2, k, nsub, index, numblk, NZmax, r, kstart, kend;

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
     numblk  = mxGetN(blk_cell_pr);            
     blksize = mxGetPr(blk_cell_pr);
     cumblksize = mxCalloc(numblk+1,sizeof(int)); 
     blknnz = mxCalloc(numblk+1,sizeof(int)); 
     cumblksize[0] = 0; blknnz[0] = 0; 
     n = 0;  n2 = 0; 
     for (k=0; k<numblk; ++k) {
          nsub = (int) blksize[k];
          n  += nsub; 
          n2 += nsub*nsub;  
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
         isspA  = mxIsSparse(A_cell_pr);
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
         if (isspA) { irA = mxGetIr(prhs[1]); 
                      jcA = mxGetJc(prhs[1]); }
         if (iscmpA) { AI = mxGetPi(prhs[1]); }
     }
     if (numblk > 1) 
        { isspB = 1; }
     else { 
        if (nrhs > 2) {isspB = (int)*mxGetPr(prhs[2]);} else {isspB = isspA;} 
     }
     if (nrhs > 4) {colidx = (int)*mxGetPr(prhs[4]) -1;} else {colidx = 0;} 
     if (colidx > n1) { 
         mexErrMsgTxt("mexsmat: colidx exceeds size(Avec,2)."); 
     }    
     /***** create return argument *****/
     if (isspB) {
        if (numblk == 1 & isspA) {           
	  NZmax = jcA[colidx+1]-jcA[colidx]; 
        } else {
	  NZmax = blknnz[numblk]; 
        }
        if (iscmpA) { 
           plhs[0] = mxCreateSparse(n,n,NZmax,mxCOMPLEX); 
        } else {
           plhs[0] = mxCreateSparse(n,n,NZmax,mxREAL); 
	}
        B = mxGetPr(plhs[0]);        
        irB = mxGetIr(plhs[0]); 
        jcB = mxGetJc(plhs[0]); 
        if (iscmpA) { BI = mxGetPi(plhs[0]); }
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
     if (numblk == 1) { 
        if (iscmpA) { 
           mat1cmp(n,A,irA,jcA,isspA,m1,colidx,B,irB,jcB,isspB,AI,BI);          
	} else {
           mat1(n,A,irA,jcA,isspA,m1,colidx,B,irB,jcB,isspB);  
	}
     } else {
       if (isspA) {
          Atmp = mxCalloc(blknnz[numblk],sizeof(double)); 
          kstart = jcA[colidx]; kend = jcA[colidx+1]; 
          for (k=kstart; k<kend; k++) { r = irA[k]; Atmp[r] = A[k]; } 
	  if (iscmpA) {  
 	     AItmp = mxCalloc(blknnz[numblk],sizeof(double)); 
             for (k=kstart; k<kend; k++) { r = irA[k]; AItmp[r] = AI[k]; } 
             mat2cmp(n,numblk,cumblksize,blknnz,Atmp,m1,0,B,irB,jcB,isspB,AItmp,BI);
	  } else { 
             mat2(n,numblk,cumblksize,blknnz,Atmp,m1,0,B,irB,jcB,isspB); 
	  }
       } else {
          if (iscmpA) { 
             mat2cmp(n,numblk,cumblksize,blknnz,A,m1,colidx,B,irB,jcB,isspB,AI,BI); 
	  } else {
             mat2(n,numblk,cumblksize,blknnz,A,m1,colidx,B,irB,jcB,isspB); 
	  }
       }
     }
     if ((isspA) & (numblk>1)) { 
        mxFree(Atmp); 
        if (iscmpA) { mxFree(AItmp); }
     }
 return;
 }
/**********************************************************/


