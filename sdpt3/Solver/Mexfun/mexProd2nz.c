/******************************************************************
* mexProd2nz: compute the elements of P whose positions are
*             specified by the variable list.    
* 
*   P = mexProd2nz(blk,A,B,list), where P = A'*B. 
*   
*   output: P is sparse, where 
*           P(i,j) = <ith column of A, jth column of B>
*
*   The variable list must have its 2nd column sorted in ascending
*   order. The first and second column of list indicate the row 
*   and column indices of elements of P to be computed. 
*    
*   E.g. list = [1 1]
*               [2 2]
*               [3 2]
*               [2 3]
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
* realdot: 
**********************************************************/
double realdot(const double *x, const int idx1, const double *y, 
               const int idx2,  const int n)

{  int i;
   double r; 

   r=0.0;
   for(i=0; i< n-3; i++){                 /* LEVEL 4 */
      r+= x[i+idx1] * y[i+idx2]; i++;
      r+= x[i+idx1] * y[i+idx2]; i++;
      r+= x[i+idx1] * y[i+idx2]; i++;
      r+= x[i+idx1] * y[i+idx2]; 
   }
   if(i < n-1){                           /* LEVEL 2 */
      r+= x[i+idx1] * y[i+idx2]; i++;
      r+= x[i+idx1] * y[i+idx2]; i++;
   }
   if(i < n){                             /* LEVEL 1 */
      r+= x[i+idx1] * y[i+idx2]; 
   }
   return r; 
}
/**********************************************************
*  B dense
*  A = nxp, B = nxp
*  P(i,j) = <ith column of A, jth column of B>
**********************************************************/
void  prod1(int m, int n, int p, 
            double *A,  mwIndex *irA, mwIndex *jcA, int isspA, 
            double *B,  double *P,  mwIndex *irP, mwIndex *jcP, 
            int *list1, int *list2, int len)

{ int  j, k, r, t, rn, jn, jold, kstart, kend, idx, count;
  double  tmp;  

     jold = -1; count = 0; 
     for (t=0; t<len; ++t) {
         r = list1[t];
         j = list2[t];  
         if (j != jold) { jn = j*n; jold = j; }       
         if (!isspA) {
            rn = r*n;
            tmp = realdot(A,rn,B,jn,n);  }
         else {
            tmp = 0;  
            kstart = jcA[r];   
            kend   = jcA[r+1]; 
            for (k=kstart; k<kend; ++k) {
                idx = irA[k]; 
                tmp += A[k]*B[idx+jn]; } 
         }
         P[count] = tmp; 
         irP[count] = r;  jcP[j+1]++;  count++; 
     }
     for (k=0; k<p; k++) { jcP[k+1] += jcP[k]; } 
     return;
}
/**********************************************************
*  B sparse.
*  A = nxm, B = nxp  
*  P(i,j) = <ith column of A, jth column of B>
**********************************************************/
void  prod2(int m, int n, int p,  
            double *A, mwIndex *irA, mwIndex *jcA, int isspA,
            double *B, mwIndex *irB, mwIndex *jcB, int isspB,
            double *P, mwIndex *irP, mwIndex *jcP, double *Btmp, 
            int *list1, int *list2, int len)  

{ int  j, k, r, t, rn, jn, jold, kstart, kend, idx, count;
  double  tmp;  

     jold = -1; count = 0;  
     for (t=0; t<len; ++t) {
         r = list1[t];  
         j = list2[t];  
         if (j != jold) {
            jn = j*n;  
            /***** copy j-th column of sparse B to Btmp *****/ 
            for (k=0; k<n; ++k) { Btmp[k] = 0; }
            kstart = jcB[j];  kend = jcB[j+1];
            for (k=kstart; k<kend; ++k) { idx = irB[k]; Btmp[idx] = B[k];  }
            jold = j; 
         }
         if (!isspA) {
            rn = r*n; 
            tmp = realdot(A,rn,Btmp,0,n); }
         else {
            tmp = 0;  
            kstart = jcA[r];   kend = jcA[r+1]; 
            for (k=kstart; k<kend; ++k) {
                idx = irA[k]; 
                tmp += A[k]*Btmp[idx]; } 
         }
         P[count] = tmp; 
         irP[count] = r;  jcP[j+1]++;   count++; 
     }
     for (k=0; k<p; k++) { jcP[k+1] += jcP[k]; } 
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
   int      *list1, *list2;
   double   *Btmp, *listtmp; 
   int       m1, n1, m2, n2, mlist, nlist, isspA, isspB, k;

/* Check for proper number of arguments */
   if (nrhs < 4){
      mexErrMsgTxt("mexProd2nz: requires at least 4 input arguments."); }
   else if (nlhs>2){ 
      mexErrMsgTxt("mexProd2nz: requires 1 output argument."); }
   
/***********************************************/     

       A = mxGetPr(prhs[1]); 
       isspA = mxIsSparse(prhs[1]);
       m1 = mxGetM(prhs[1]); 
       n1 = mxGetN(prhs[1]); 
       if (isspA) { irA = mxGetIr(prhs[1]);
                    jcA = mxGetJc(prhs[1]); }
       B = mxGetPr(prhs[2]); 
       isspB  = mxIsSparse(prhs[2]);
       m2 = mxGetM(prhs[2]); 
       n2 = mxGetN(prhs[2]); 
       if (isspB) { irB = mxGetIr(prhs[2]);
                    jcB = mxGetJc(prhs[2]); }
       if (m1!=m2) { 
          mexErrMsgTxt("mexProd2nz: A, B are not compatible"); } 
       mlist = mxGetM(prhs[3]); 
       nlist = mxGetN(prhs[3]); 
       if (nlist != 2 ) { 
          mexErrMsgTxt("mexProd2nz: list must have 2 columns");}
       listtmp = mxGetPr(prhs[3]); 
       list1 = mxCalloc(mlist,sizeof(int)); 
       list2 = mxCalloc(mlist,sizeof(int));
       for (k=0; k<mlist; k++) {
          /** subtract 1 to adjust for Matlab index **/
          list1[k] = (int)listtmp[k] -1;
          list2[k] = (int)listtmp[mlist+k] -1;
       }
       plhs[0] = mxCreateSparse(n1,n2,mlist,mxREAL); 
       P=mxGetPr(plhs[0]); irP=mxGetIr(plhs[0]); jcP=mxGetJc(plhs[0]);
       /**********************************************
        * Do the actual computations in a subroutine 
        **********************************************/
        if (!isspB) { 
   	   prod1(n1,m2,n2,A,irA,jcA,isspA,B,P,irP,jcP,list1,list2,mlist); }
        else {
           Btmp = mxCalloc(m2,sizeof(double)); 
   	   prod2(n1,m2,n2,A,irA,jcA,isspA,B,irB,jcB,isspB,\
                 P,irP,jcP,Btmp,list1,list2,mlist);
	}
        mxFree(list1); mxFree(list2); 
        if (isspB) { mxFree(Btmp); }
return;
}		 
/**********************************************************/

