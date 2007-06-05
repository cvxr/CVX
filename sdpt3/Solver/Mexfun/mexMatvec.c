/***********************************************************************
* mexMatvec: compute
* 
* mexMatvec(A,y,options)
*          
*  options = 0, compute A*y
*          = 1, compute (y'*A)' 
*
* Copyright (c) 2004 by
* K.C. Toh
* Last Modified: 120404   
************************************************************************/

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

/********************************
* realdotde: x dense matrix,  
*            y dense vector
*********************************/
double realdotde(const double *x, const int idx, 
                 const double *y, const int n)

{  int i;
   double r; 

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
/********************************
* saxpymat:  z = z + alpha*y 
* y dense matrix, z dense vector
********************************/
void saxpymat(const double *y,  const int idx1, 
              const int istart, const int iend,
              const double alp, double *z, const int idx2)
{  int i;
  
   for(i=istart; i< iend-3; i++){             /* LEVEL 4 */
      z[i+idx2] += alp * y[i+idx1]; i++;
      z[i+idx2] += alp * y[i+idx1]; i++;
      z[i+idx2] += alp * y[i+idx1]; i++;
      z[i+idx2] += alp * y[i+idx1]; 
   }
   if(i < iend-1){                            /* LEVEL 2 */
      z[i+idx2] += alp * y[i+idx1]; i++;
      z[i+idx2] += alp * y[i+idx1]; i++;
   }
   if(i < iend){                              /* LEVEL 1 */
      z[i+idx2] += alp * y[i+idx1]; 
   }
   return;
}
/**********************************************************/
void mexFunction(int nlhs,   mxArray  *plhs[], 
                 int nrhs,   const mxArray  *prhs[] )

{    double   *A, *y;
     mwIndex  *irA, *jcA, *iry, *jcy; 
     double   *ytmp, *Ay;
     int      isspA, isspy, options;

     int      m1, n1, m2, n2, j, jm1; 
     int      i, r, k, istart, iend, kstart, kend;
     double   tmp; 

/* CHECK THE DIMENSIONS */

   if (mxIsCell(prhs[0]) | mxIsCell(prhs[1])) { 
       mexErrMsgTxt(" mexMatvec: A, x must be a double array"); }
   if (nrhs <2) {
       mexErrMsgTxt(" mexMatvec: must have at least 2 inputs"); }
   if (nlhs > 2) { 
       mexErrMsgTxt("mexMatvec: requires 1 output argument"); }
   if (nrhs == 2) { options = 0; } 
   else { options = (int) *mxGetPr(prhs[2]); } 


   /***** assign pointers *****/

       A = mxGetPr(prhs[0]); 
       m1 = mxGetM(prhs[0]); 
       n1 = mxGetN(prhs[0]);
       isspA = mxIsSparse(prhs[0]);
       if (isspA) { irA = mxGetIr(prhs[0]);
	            jcA = mxGetJc(prhs[0]); }
       isspy = mxIsSparse(prhs[1]);
       m2 = mxGetM(prhs[1]); 
       n2 = mxGetN(prhs[1]);            
       if (n2 > 1) { 
 	  mexErrMsgTxt("mexMatvec: 2ND input must be a column vector"); }
       if (isspy) { 
          iry = mxGetIr(prhs[1]);
          jcy = mxGetJc(prhs[1]); 
          ytmp = mxGetPr(prhs[1]); 
          /***** copy ytmp to y *****/ 
          y = mxCalloc(m2,sizeof(double)); 
          kstart = jcy[0]; kend = jcy[1]; 
          for (k=kstart; k<kend; k++) { 
	      r = iry[k]; y[r] = ytmp[k]; }
       } else {
          y = mxGetPr(prhs[1]); 
       }       
       if (options == 0 & n1 != m2) {
          mexErrMsgTxt("mexMatvec: 1ST and 2ND input not compatible."); 
       } else if (options & m1 != m2) {
          mexErrMsgTxt("mexMatvec: 1ST and 2ND input not compatible."); 
       }

       /***** create return argument *****/
       if (options==0) { 
          plhs[0] = mxCreateDoubleMatrix(m1,1,mxREAL); 
       } else { 
          plhs[0] = mxCreateDoubleMatrix(n1,1,mxREAL); 
       }
       Ay = mxGetPr(plhs[0]); 

    /***** main body *****/
    if (options==0) {
       if (!isspA) {
          for (j=0; j<n1; j++){ 
	     jm1 = j*m1;
             tmp = y[j]; 
             if (tmp !=0) { 
                saxpymat(A,jm1,0,m1,tmp,Ay,0); } 
	  }
       } else {
          for (j=0; j<n1; j++){
             tmp = y[j];
             if (tmp != 0) {
                istart = jcA[j]; iend = jcA[j+1]; 
  	        for (i=istart; i<iend; i++) {
                    r = irA[i]; 
	            Ay[r] += tmp*A[i]; }
	     }
	  }
       }
    } else {
       if (!isspA) {
          for (j=0; j<n1; j++){ 
            jm1 = j*m1; 
	    Ay[j] = realdotde(A,jm1,y,m1); 
	  }
       } else {
          for (j=0; j<n1; j++){
               istart = jcA[j]; iend = jcA[j+1]; 
               tmp = 0; 
  	       for (i=istart; i<iend; i++) {
                   r = irA[i]; 
	           tmp += y[r]*A[i]; }
               Ay[j] = tmp; 
	  }	  
       }
    }
 return;
}
/**********************************************************/

