/**********************************************************************
*   y = mexbsolve(U,b,options)
*  Given upper-triangular matrix U, 
*  options = 1 (default), solves     U *y = b, 
*          = 2          , solves     U'*y = b. 
*  Using loop unrolling, it is about 10 times faster than MATLAB's
*  built in solver, on Intel Pentium PRO 200 (this typically depends
*  on the system's architecture).
*
*  Note: This program is a modification of a program that is 
*        part of SeDuMi 1.02 (03AUG1998) written by Jos F. Sturm. 
*********************************************************************/

#include "mex.h"
#include <math.h>

#if !defined(SQR)
#define SQR(x) ((x)*(x))
#endif

#if !defined(MIN)
#define  MIN(A, B)   ((A) < (B) ? (A) : (B))
#endif

/**************************************************************
   TIME-CRITICAL PROCEDURE -- r=realdot(x,y,n)
   Computes r=sum(x_i * y_i) using LEVEL 8 loop-unrolling.
**************************************************************/
 double realdot(const double *x, const double *y, const int n)
 {
   int i;
   double r;

   r=0.0;
   for(r=0.0, i=0; i< n-7; i++){          /* LEVEL 8 */
      r+= x[i] * y[i]; i++;
      r+= x[i] * y[i]; i++;
      r+= x[i] * y[i]; i++;
      r+= x[i] * y[i]; i++;
      r+= x[i] * y[i]; i++;
      r+= x[i] * y[i]; i++;
      r+= x[i] * y[i]; i++;
      r+= x[i] * y[i];
   }
   if(i < n-3){                           /* LEVEL 4 */
     r+= x[i] * y[i]; i++;
     r+= x[i] * y[i]; i++;
     r+= x[i] * y[i]; i++;
     r+= x[i] * y[i]; i++;
   }
   if(i < n-1){                           /* LEVEL 2 */
     r+= x[i] * y[i]; i++;
     r+= x[i] * y[i]; i++;
   }
   if(i < n)                              /* LEVEL 1 */
     r+= x[i] * y[i];
   return r;
}
/*************************************************************
   TIME-CRITICAL PROCEDURE -- subscalarmul(x,alpha,y,n)
   Computes x -= alpha * y using LEVEL 8 loop-unrolling.
**************************************************************/
void subscalarmul(double *x, const double alpha, 
                  const double *y, const int n)
{
 int i;

  for(i=0; i< n-7; i++){          /* LEVEL 8 */
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i];
  }
  if(i < n-3){                    /* LEVEL 4 */
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
  }
  if(i < n-1){                    /* LEVEL 2 */
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
  }
  if(i < n)                       /* LEVEL 1 */
    x[i] -= alpha * y[i];
}
/*************************************************************
   PROCEDURE lbsolve -- Solve y from U'*y = x. 
   INPUT
     u - n x n full matrix with only upper-triangular entries
     x,n - length n vector.
   OUTPUT
     y - length n vector, y = U'\x.
**************************************************************/
void lbsolve(double *y,const double *u,const double *x,const int n)
{
  int k;

  /*------------------------------------------------------------
     The first equation, u(:,1)'*y=x(1), yields y(1) = x(1)/u(1,1).
     For k=2:n, we solve
        u(1:k-1,k)'*y(1:k-1) + u(k,k)*y(k) = x(k).
   -------------------------------------------------------------*/
  for (k = 0; k < n; k++, u += n)
      y[k] = (x[k] - realdot(y,u,k)) / u[k];
}
/*************************************************************
   PROCEDURE ubsolve -- Solves xnew from U * xnew = x,
     where U is upper-triangular.
   INPUT
     u,n - n x n full matrix with only upper-triangular entries
   UPDATED
     x - length n vector
     On input, contains the right-hand-side.
     On output, xnew = U\xold
**************************************************************/
void ubsolve(double *x,const double *u,const int n)
{
   int j;

  /*------------------------------------------------------------
     At each step j= n-1,...0, we have a (j+1) x (j+1) upper triangular
     system "U*xnew = x". The last equation gives:
       xnew[j] = x[j] / u[j,j]
     Then, we update the right-hand side:
       xnew(0:j-1) = x(0:j-1) - xnew[j] * u(0:j-1)
   --------------------------------------------------------------*/
   j = n;
   u += SQR(n);
   while(j > 0){
     --j;
     u -= n;
     x[j] /= u[j];
     subscalarmul(x,x[j],u,j);
   }
}
/*************************************************************
*   PROCEDURE mexFunction - Entry for Matlab
**************************************************************/
 void mexFunction(const int nlhs, mxArray *plhs[],
                  const int nrhs, const mxArray *prhs[])
{
   const double *u, *b;
   const int    *irb, *jcb; 
   int    n, k, kend, isspb, options;   
   double *y, *btmp;

   if (nrhs < 2) {
      mexErrMsgTxt("mexbsolve requires 2 input arguments."); }
   if (nlhs > 1) {
      mexErrMsgTxt("mexbsolve generates 1 output argument."); }
 
   u = mxGetPr(prhs[0]);
   if (mxIsSparse(prhs[0])) {
      mexErrMsgTxt("Sparse U not supported by mexbsolve."); }
   n = mxGetM(prhs[0]); 
   if (mxGetN(prhs[0]) != n) {
      mexErrMsgTxt("U should be square and upper triangular."); }
   b = mxGetPr(prhs[1]);
   isspb = mxIsSparse(prhs[1]); 
   if ( mxGetM(prhs[1])*mxGetN(prhs[1]) != n ) {
       mexErrMsgTxt("Size mismatch in U,b."); }
   if (nrhs > 2) { 
      options = (int)*mxGetPr(prhs[2]); }
   else {
      options = 1; }

   plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
   y = mxGetPr(plhs[0]);
  
   if (isspb) {
      btmp = mxCalloc(n,sizeof(double));       
      irb = mxGetIr(prhs[1]); jcb = mxGetJc(prhs[1]); 
      kend = jcb[1]; 
      for (k=0; k<kend; k++) { btmp[irb[k]] = b[k]; } 
   }
   if (options==1) { 
      if (isspb) { 
         for (k=0; k<n; k++) { y[k]=btmp[k]; } }
      else {
         for (k=0; k<n; k++) { y[k]=b[k]; } }
      ubsolve(y,u,n);  
   }
   else if (options==2) { 
       if (isspb) { lbsolve(y,u,btmp,n); }
       else       { lbsolve(y,u,b,n); }
   }   
   if (isspb) { mxFree(btmp); }
}
/*************************************************************/
