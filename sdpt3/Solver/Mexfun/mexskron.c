/***********************************************************************
* mexskron.c : *  Find the matrix presentation of 
*  symmetric kronecker product skron(A,B), where 
*  A,B are symmetric. 
*
*  K = mexskron(blk,A,B,sym); 
*
*  sym = 1 if A=B
*      = 0 otherwise
*
* (ij)-column = 0.5*svec(AUB + BUA)*xij,  where 
*       xij = 1/2       if i=j 
*           = 1/sqrt(2) otherwise
*       U   = ei*ej' + ej*ei'. 
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

#if !defined(r2)
#define  r2   1.41421356237309504880      /* sqrt(2) */
#endif

#if !defined(ir2)
#define  ir2  0.70710678118654752440      /* 1/sqrt(2) */
#endif

/**********************************************************
* A, B may not be equal
***********************************************************/
void skron(int n, int maxblksize, double *P, double *Q,
           double *x1, double *y1, double *x2, double *y2,
           int r, int c, double *vvtmp)

{  int idx, i, j, k, rn, cn;
   double tmp, tmp1, tmp2, tmp3, tmp4; 
   double hf=0.5;

   rn = r*maxblksize; cn = c*maxblksize;
   for (k=0; k<n; k++) {
      x1[k] = P[k+rn]; y1[k] = Q[k+cn];  
      x2[k] = P[k+cn]; y2[k] = Q[k+rn];
   } 
   if (r < c) {
      idx = 0;
      for (j=0; j<n; j++) {
          tmp1 = hf*y1[j]; tmp2 = hf*y2[j]; 
          tmp3 = hf*x1[j]; tmp4 = hf*x2[j];   
          for (i=0; i<j; i++) {
              vvtmp[idx] = tmp1*x1[i]+tmp2*x2[i]+tmp3*y1[i]+tmp4*y2[i];
              idx++; }
          tmp = tmp1*x1[j]+tmp2*x2[j]+tmp3*y1[j]+tmp4*y2[j];
          vvtmp[idx] = tmp*ir2; 
          idx++;
      }
   } else {   
      /*** r = c ***/
      idx = 0;
      for (j=0; j<n; j++) {
  	  tmp1 = ir2*y1[j]; tmp2 = ir2*x1[j]; 
          for (i=0; i<j; i++) {
              vvtmp[idx] = tmp1*x1[i]+tmp2*y1[i];
              idx++; }
          vvtmp[idx] = y1[j]*x1[j];
          idx++;
      }
   }
return;
}
/**********************************************************
* A=B 
***********************************************************/
void skron2(int n, int maxblksize, double *P, double *Q,
           double *x1, double *y1, double *x2, double *y2,
           int r, int c, double *vvtmp)

{  int idx, i, j, k, rn, cn;
   double tmp1, tmp2, tmp; 
 
   rn = r*maxblksize; cn = c*maxblksize;
   for (k=0; k<n; k++) {
      x1[k] = P[k+rn]; y1[k] = Q[k+cn];  
      x2[k] = P[k+cn]; y2[k] = Q[k+rn];
   } 
   if (r < c) {
      idx = 0;
      for (j=0; j<n; j++) {
          tmp1 = y1[j]; tmp2 = y2[j]; 
          for (i=0; i<j; i++) {
              vvtmp[idx] = tmp1*x1[i] + tmp2*x2[i];
              idx++; }
          tmp = tmp1*x1[j] + tmp2*x2[j];
          vvtmp[idx] = tmp*ir2; 
          idx++;
      }
   } else {   
      /*** r = c ***/
      idx = 0;
      for (j=0; j<n; j++) {
          tmp1 = r2*y1[j]; 
          for (i=0; i<j; i++) {
              vvtmp[idx] = tmp1*x1[i];
              idx++; }
          vvtmp[idx] = y1[j]*x1[j];
          idx++;
      }
   }
return;
}
/**********************************************************
* 
***********************************************************/
void mexFunction(int nlhs, mxArray  *plhs[], 
                 int nrhs, const mxArray  *prhs[] )

{    mxArray  *blk_cell_pr;
     mxArray  *tmparr[3];
     double   *A,  *B,  *blksizetmp, *P, *Q, *V; 
     double   *ii, *jj, *vv, *vvtmp, *x1, *y1, *x2, *y2; 
     mwIndex  *irA, *jcA, *irB, *jcB, *irV, *jcV;
     int      *blksize, *blksize2, *cumblksize, *blksize4;
     int       isspA, isspB, sym;

     mwIndex   subs[2];
     mwSize    nsubs=2; 
     int    m, n, n2, nsub, k, index, numblk, len, maxblksize; 
     int    i, j, l, jn, idxstart, idxend, kstart, kend, idx, blklen; 
     int    blklen2, rowidx, colidx; 

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs < 3){
      mexErrMsgTxt("mexskron: requires at least 3 input arguments."); }
   if (nlhs != 1){ 
      mexErrMsgTxt("mexskron: requires 1 output argument."); }

/* CHECK THE DIMENSIONS */

    if (mxGetM(prhs[0]) > 1) { 
       mexErrMsgTxt("mexskron: blk can only have 1 row."); }
    if (mxGetN(prhs[0]) != 2) { 
       mexErrMsgTxt("mexskron: blk must have 2 columns."); }
    subs[0] = 0; subs[1] = 1;
    index = mxCalcSingleSubscript(prhs[0],nsubs,subs); 
    blk_cell_pr = mxGetCell(prhs[0],index);
    numblk  = mxGetN(blk_cell_pr);
    blksizetmp = mxGetPr(blk_cell_pr); 
    blksize    = mxCalloc(numblk,sizeof(int));
    for (k=0; k<numblk; k++) { blksize[k] = (int) blksizetmp[k]; }
    cumblksize = mxCalloc(numblk+1,sizeof(int));
    blksize2   = mxCalloc(numblk+1,sizeof(int)); 
    blksize4   = mxCalloc(numblk+1,sizeof(int)); 
    cumblksize[0] = 0; blksize2[0] = 0;  blksize4[0] = 0; 
    maxblksize=0; 
    for (k=0; k<numblk; ++k) {
        nsub = blksize[k];
        n2   = nsub*(nsub+1)/2;
        cumblksize[k+1] = cumblksize[k] + nsub;                 
        blksize2[k+1]   = blksize2[k] + n2; 
        blksize4[k+1]   = blksize4[k] + n2*n2; 
        maxblksize = MAX(maxblksize,nsub);  
    }
    /***** assign pointers *****/
    m = mxGetM(prhs[1]); 
    n = mxGetN(prhs[1]); 
    if (m != n) { 
       mexErrMsgTxt("mexskron: matrix A must be square."); }
    A = mxGetPr(prhs[1]); 
    isspA = mxIsSparse(prhs[1]); 
    if (isspA) { irA = mxGetIr(prhs[1]); 
                 jcA = mxGetJc(prhs[1]); }             
    m = mxGetM(prhs[2]); 
    n = mxGetN(prhs[2]); 
    if (m != n) { 
       mexErrMsgTxt("mexskron: matrix B must be square."); }
    B = mxGetPr(prhs[2]); 
    isspB = mxIsSparse(prhs[2]); 
    if (isspB) { irB = mxGetIr(prhs[2]); 
                 jcB = mxGetJc(prhs[2]); }             
    if ((isspA+isspB)==1) {
       mexErrMsgTxt("mexskron: A and B must be both sparse or both dense"); 
    }
    sym = 0; 
    if (nrhs == 4) {
       if (mxGetPr(prhs[3]) != NULL) { 
          sym = (int)*mxGetPr(prhs[3]); }
    }
    /***** create temporary array *****/
    tmparr[0] = mxCreateDoubleMatrix(maxblksize,maxblksize,mxREAL);
    P = mxGetPr(tmparr[0]);
    tmparr[1] = mxCreateDoubleMatrix(maxblksize,maxblksize,mxREAL);
    Q = mxGetPr(tmparr[1]);
    tmparr[2] = mxCreateDoubleMatrix(maxblksize*(maxblksize+1)/2,1,mxREAL); 
    vvtmp = mxGetPr(tmparr[2]);

    x1 = mxCalloc(maxblksize,sizeof(double)); 
    y1 = mxCalloc(maxblksize,sizeof(double)); 
    x2 = mxCalloc(maxblksize,sizeof(double)); 
    y2 = mxCalloc(maxblksize,sizeof(double)); 

    /***** create return argument *****/
    len = blksize4[numblk]; 
    plhs[0] = mxCreateSparse(blksize2[numblk],blksize2[numblk],len,mxREAL); 
    V = mxGetPr(plhs[0]);  
    irV = mxGetIr(plhs[0]); 
    jcV = mxGetJc(plhs[0]); 

    /***** Do the computations in a subroutine *****/

    for (l=0; l<numblk; l++) {
        blklen = blksize[l];
        /*----- copy lth block to P,Q -----*/ 
        for (j=0; j<blklen; j++) {
	    jn = j*maxblksize; 
	    for (k=0; k<blklen; k++) {
	        idx = k+jn; 
                P[idx] = 0; Q[idx] = 0; }	    
        }
        idxstart = cumblksize[l]; 
        idxend   = cumblksize[l+1]; 
        if (isspA & isspB) {
           for (j=idxstart; j<idxend; j++) {
	       kstart = jcA[j]; kend = jcA[j+1];
               jn = (j-idxstart)*maxblksize; 
               for (k=kstart; k<kend; k++) {
	           idx = irA[k]-idxstart; 
	           P[idx+jn] = A[k]; }	
               kstart = jcB[j]; kend = jcB[j+1]; 
               for (k=kstart; k<kend; k++) {
	           idx = irB[k]-idxstart; 
	           Q[idx+jn] = B[k]; }	   
	   }
	} else {
           for (j=idxstart; j<idxend; j++) {
	       colidx = j*cumblksize[numblk]; 
               jn = (j-idxstart)*maxblksize; 
               for (k=idxstart; k<idxend; k++) {
		   idx = (k-idxstart) + jn;
		   P[idx] = A[k+colidx]; 
	           Q[idx] = B[k+colidx]; 	 
               }
	   }
	}
        /*----- skron(P,Q) -----*/ 
        blklen2 = blklen*(blklen+1)/2; 
        idx = blksize4[l];
        for (j=0; j<blklen; j++) {
	   for (i=0; i<=j; i++) {
	      if (sym==0) {
	        skron(blklen,maxblksize,P,Q,x1,y1,x2,y2,i,j,vvtmp);
	      } else {
	        skron2(blklen,maxblksize,P,Q,x1,y1,x2,y2,i,j,vvtmp);
	      }
              colidx = blksize2[l] + i + j*(j+1)/2; 
              rowidx = blksize2[l];
	      for (k=0; k<blklen2; k++) {            
		 V[idx] = vvtmp[k];
                 irV[idx] = rowidx+k; 
                 idx++;
	      }
              jcV[colidx+1] = idx;         
           }
	}  
    }
    /***** free memory *************************/  
    mxFree(blksize);    mxFree(blksize2); 
    mxFree(cumblksize); mxFree(blksize4); 
    mxFree(x1); mxFree(y1); mxFree(x2); mxFree(y2); 
    mxDestroyArray(*tmparr);   
return;
}
/**********************************************************/

