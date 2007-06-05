/***********************************************************************
* mexqops.c : C mex file 
*
*   z = mexqops(blk,x,y,options);
*
*   Input: blk   = [n1, n2, ... nk]
*          x = n-vector or k-vector
*          y = n-vector, where n = n1+...+nk
*  
*   options = 1, z = k-vector, z(i) = <xi,yi> 
*           = 2, z = k-vector, z(i) = 2*xi(1)yi(1) - <xi,yi>
*           = 3, z = n-vector, zi   = x(i)*yi 
*           = 4, z = n-vector, zi   = x(i)*yi, zi(1) = -zi(1). 
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
* ops1 
**********************************************************/
void ops1(double *x, double *y, double *z, 
          int numblk, int *cumblk, int options) 

{  int  j, l, jstart, jend;
   double tmp; 
   
   if (options == 1) { 
      for (l=0; l<numblk; ++l) {  
          jstart = cumblk[l]; jend = cumblk[l+1];
          tmp = 0; 
          for (j=jstart; j<jend; ++j) { tmp += x[j]*y[j]; }
          z[l] = tmp; 
      }
   }
   else if (options == 2) { 
      for (l=0; l<numblk; ++l) {  
          jstart = cumblk[l]; jend = cumblk[l+1];
          tmp = x[jstart]*y[jstart]; 
          for (j=jstart+1; j<jend; ++j) { tmp -= x[j]*y[j]; }
          z[l] = tmp;
      }          
   }
   return;
}
/**********************************************************
* ops3 
**********************************************************/
void ops3(double *x, double *y, double *z, 
          int numblk, int *cumblk, int options) 

{  int  j, l, jstart, jend;
   double tmp;
   
   if (options == 3) {
      for (l=0; l<numblk; ++l) {  
         jstart = cumblk[l]; jend = cumblk[l+1];
         tmp = x[l];
         for (j=jstart; j<jend; ++j) { 
             z[j] = tmp*y[j]; }
      }
   }
   else if (options == 4) { 
      for (l=0; l<numblk; ++l) {  
         jstart = cumblk[l]; jend = cumblk[l+1];
         tmp = x[l];
         for (j=jstart; j<jend; ++j) { 
             z[j] = tmp*y[j]; }
         z[jstart] = -z[jstart];
      }
   }
   return;
}
/**********************************************************
* 
***********************************************************/
void mexFunction(int nlhs, mxArray  *plhs[], 
                 int nrhs, const mxArray  *prhs[] )

{    double   *blksize, *x, *y, *z, *xtmp, *ytmp; 
     mwIndex  *irx, *jcx, *iry, *jcy; 
     int      *cumblk;
     int   mblk, options; 

     int   l, n, k, r, numblk, cols;

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs !=4){
      mexErrMsgTxt("mexqops: requires 4 input arguments."); }
   if (nlhs > 1){ 
      mexErrMsgTxt("mexqops: requires 1 output argument."); }

/* CHECK THE DIMENSIONS */
    
    numblk = mxGetN(prhs[0]); 
    blksize = mxGetPr(prhs[0]);
    cumblk = mxCalloc(numblk+1,sizeof(int)); 
    for (l=0; l<numblk; l++) { 
        cols = (int)blksize[l]; 
        cumblk[l+1] = cumblk[l] + cols;  
    } 
    n = cumblk[numblk]; 
    /***** assign pointers *****/
    options = (int)*mxGetPr(prhs[3]); 
    if (mxGetM(prhs[2]) != n) {
       mexErrMsgTxt("mexqops: dim not compatible.");  }
    if (options < 3 & mxGetM(prhs[1]) != n) {        
       mexErrMsgTxt("mexqops: dim not compatible.."); }            
    if (options >= 3 & mxGetM(prhs[1]) != numblk) { 
       mexErrMsgTxt("mexqops: dim not compatible..."); }    
    if (mxIsSparse(prhs[1])) { 
       irx = mxGetIr(prhs[1]);  jcx = mxGetJc(prhs[1]); xtmp = mxGetPr(prhs[1]); 
       x = mxCalloc(n,sizeof(double)); 
       for (k=0; k<jcx[1]; ++k) { r=irx[k]; x[r]=xtmp[k]; } }
    else { 
       x = mxGetPr(prhs[1]);    
    }    
    if (mxIsSparse(prhs[2])) { 
       iry = mxGetIr(prhs[2]);  jcy = mxGetJc(prhs[2]); ytmp = mxGetPr(prhs[2]); 
       y = mxCalloc(n,sizeof(double)); 
       for (k=0; k<jcy[1]; ++k) { r=iry[k]; y[r]=ytmp[k]; }  }
    else { 
       y = mxGetPr(prhs[2]); 
    }
    /***** Do the computations in a subroutine *****/
    if (options < 3) {        
       plhs[0] = mxCreateDoubleMatrix(numblk,1,mxREAL);
       z = mxGetPr(plhs[0]);
       ops1(x,y,z,numblk,cumblk,options); 
    }
    else if (options >= 3) { 
       plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
       z = mxGetPr(plhs[0]);
       ops3(x,y,z,numblk,cumblk,options); }
 return;
 }
/**********************************************************/
