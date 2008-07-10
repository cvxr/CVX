/*
%                                 ddotX = ddot(d,X,blkstart [, Xblkjc])
% DDOT Given N x m matrix X, creates (blkstart(end)-blkstart(1)) x m matrix
%   ddotX, having entries d[i]'* xj[i] for each (Lorentz norm bound) block
%   blkstart(i):blkstart(i+1)-1. If X is sparse, then Xblkjc(:,2:3) should
%   point to first and 1-beyond-last nonzero in blkstart range for each column.
%
% SEE ALSO sedumi, partitA.
% **********  INTERNAL FUNCTION OF SEDUMI **********
function ddotX = ddot(d,X,blkstart, Xblkjc)

% This file is part of SeDuMi 1.1 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2005 McMaster University, Hamilton, CANADA  (since 1.1)
%
% Copyright (C) 2001 Jos F. Sturm (up to 1.05R5)
%   Dept. Econometrics & O.R., Tilburg University, the Netherlands.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% Affiliation SeDuMi 1.03 and 1.04Beta (2000):
%   Dept. Quantitative Economics, Maastricht University, the Netherlands.
%
% Affiliations up to SeDuMi 1.02 (AUG1998):
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA


*/

#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define DDOTX_OUT plhs[0]
#define NPAROUT 1

#define D_IN prhs[0]
#define X_IN prhs[1]
#define BLKSTART_IN prhs[2]
#define NPARINMIN 3
#define XBLKJC_IN prhs[3]
#define NPARIN 4

/* ************************************************************
   PROCEDURE ddotxj -Compute y[k]= d[k]'*xpr[k] for each lorentz block k.
   INPUT
     d - qDim scaling vector with qDim := blkstart[nblk]-blkstart[0].
     xpr - qDim data vector. 
     blkstart - length nblk+1 array, listing 1st subscript per block.
       NOTE: should have blkstart[0] == 0.
     nblk - Number of blocks.
   OUTPUT
     ypr - nblk vector. Gives d[k]'*xj[k] for each block.
   ************************************************************ */
void ddotxj(double *ypr, const double *d, const double *xpr,
            const mwIndex *blkstart, const mwIndex nblk)
{
  mwIndex k;
  mxAssert(blkstart[0] == 0,"");
  for(k = 0; k < nblk; k++)
    ypr[k] = realdot(d+blkstart[k],xpr+blkstart[k], blkstart[k+1]-blkstart[k]);
}

/* ************************************************************
   PROCEDURE spddotxj - Compute y[k] = d_k'*xj_k for each nonzero
      block in xj.
   INPUT
     d - qDim scaling vector with qDim := blkstart[nblk]-blkstart[0].
     xir, xpr - sparse matrix. We compute d[k]'*xj[k] for each (lorentz) block
       where the column xj has nonzeros.
     xjc0, xjc1 - Length m arrays, subscripts of column j in blkstart-range
        are between xjc0(j) and xjc1(j).
     blkstart - length nblk+1 array. Lorentz block k has subscripts
       blkstart[k]:blkstart[k+1]-1.
     xblk - length qDim array, with k = xblk(i-blkstart[0]) iff
       blkstart[k] <= i < blkstart[k+1], k=0:nblk-1.
   OUTPUT
     y - sparse nblk x m matrix, with y.jc[m] <= sum(xjc1-xjc0).
        y(k,j) = d[k]'*xj[k]
   ************************************************************ */
void spddotxj(jcir y, const double *d,
              const mwIndex *xir, const double *xpr, const mwIndex *xjc0,
              const mwIndex *xjc1, const mwIndex *xblk, const mwIndex *blkstart,
              const mwIndex nblk, const mwIndex m)
{
  mwIndex knz, nexti, inz, i, j, k, lend;
  double yk;
/* ------------------------------------------------------------
   INIT: Let blkstart[0] point to 1st nonzero in d and xblk, and
   let knz poin to 1st available entry in y.
   Let lend := blkstart[nblk] be 1 beyond valid subscripts.
   ------------------------------------------------------------ */
  d -= blkstart[0];           /* Make d=d(blkstart[0]:blkstart[lorN]) */
  xblk -= blkstart[0];
  knz = 0;
  lend = blkstart[nblk];
  for(j = 0; j < m; j++){
    y.jc[j] = knz;
/* ------------------------------------------------------------
   Process column only if nonzero subscripts in blkstart[0:nblk].
   ------------------------------------------------------------ */
    if((inz = xjc0[j]) < xjc1[j])
      if( (i = xir[inz]) < lend){
/* ------------------------------------------------------------
   Open initial block k; current block has subscripts smaller than nexti.
   Accumulate yk = ddotxj[k].
   ------------------------------------------------------------ */
        k = xblk[i];
        nexti = blkstart[k + 1];
        yk = d[i] * xpr[inz];
/* ------------------------------------------------------------
   Browse through nonzeros in xj
   ------------------------------------------------------------ */
        for(++inz; inz < xjc1[j]; inz++)
          if( (i = xir[inz]) < nexti)
            yk += d[i] * xpr[inz];
          else if(i < lend){
/* ------------------------------------------------------------
   If we finished the previous nonzero Lorentz block, then write entry,
   and initialize new block.
   ------------------------------------------------------------ */
            y.ir[knz] = k;              /* yir lists Lorentz blocks */
            y.pr[knz++] = yk;
            k = xblk[i];               /* init new Lorentz block */
            nexti = blkstart[k + 1];
            yk = d[i] * xpr[inz];
          }
          else                   /* finished with all Lorentz blocks */
            break;
/* ------------------------------------------------------------
   Write last yk = ddotxj[k] entry into y(:,j).
   ------------------------------------------------------------ */
        y.ir[knz] = k;              /* yir lists Lorentz blocks */
        y.pr[knz++] = yk;
      } /* If column j has valid nonzeros */
  } /* j=0:m-1 */
/* ------------------------------------------------------------
   Close last column of y
   ------------------------------------------------------------ */
  y.jc[m] = knz;
}

/* ============================================================
   MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  mwIndex i, j, k, m, nrows, maxnnz, nblk, qDim;
  const double *d, *XjcPr, *blkstartPr;
  mwIndex *xjc1, *xblk, *blkstart;
  jcir X, ddotx;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARINMIN, "ddot requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "ddot generates less output arguments.");
/* ------------------------------------------------------------
   Get INPUTS d, X, blkstart.
   ------------------------------------------------------------ */
  d = mxGetPr(D_IN);
  qDim = mxGetM(D_IN) * mxGetN(D_IN);
  nrows = mxGetM(X_IN);
  m = mxGetN(X_IN);
  X.pr = mxGetPr(X_IN);
  blkstartPr = mxGetPr(BLKSTART_IN);
  nblk = mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN) - 1;
  mxAssert(nblk >= 0, "blkstart size mismatch.");
/* ------------------------------------------------------------
   Allocate mwIndex working array blkstart(nblk+1).
   ------------------------------------------------------------ */
  blkstart = (mwIndex *) mxCalloc(nblk + 1, sizeof(mwIndex));
/* ------------------------------------------------------------
   Convert Fortran double to C mwIndex
   ------------------------------------------------------------ */
  for(i = 0; i <= nblk; i++){
    j = (mwIndex) blkstartPr[i];             /* double to mwIndex */
    mxAssert(j>0,"");
    blkstart[i] = --j;
  }
  if(qDim != blkstart[nblk] - blkstart[0]){
    mxAssert(qDim >= blkstart[nblk], "d size mismatch.");
    d += blkstart[0];                   /* Point to Lorentz norm-bound */
    qDim = blkstart[nblk] - blkstart[0];
  }
/* ------------------------------------------------------------
   CASE THAT X IS FULL:
   ------------------------------------------------------------ */
  if(!mxIsSparse(X_IN)){
    if(nrows != qDim)
      if(nrows < blkstart[nblk]){
         mxAssert(nrows == nblk + qDim, "X size mismatch");
         X.pr += nblk;                 /* Lorentz tr + norm bound */
      }
      else                        /* LP, Lorentz, PSD */
	X.pr += blkstart[0];      /* Point to Lorentz norm-bound */
/* ------------------------------------------------------------
   DDOTX is full nblk x m.
   ------------------------------------------------------------ */
    DDOTX_OUT = mxCreateDoubleMatrix(nblk, m, mxREAL);
    ddotx.pr = mxGetPr(DDOTX_OUT);
/* ------------------------------------------------------------
   Let blkstart -= blkstart[0], so that blkstart[0] = 0.
   ------------------------------------------------------------ */
    j = blkstart[0];
    for(i = 0; i <= nblk; i++)
      blkstart[i] -= j;
/* ------------------------------------------------------------
   Compute d[k]'*x[k,i] for all Lorentz blocks k.
   ------------------------------------------------------------ */
    for(i = 0; i < m; i++){
      ddotxj(ddotx.pr, d, X.pr, blkstart, nblk);
      ddotx.pr += nblk;
      X.pr += nrows;                /* to next column */
    }
  }
  else{
/* ------------------------------------------------------------
   The CASE that  X is SPARSE:
   ------------------------------------------------------------ */
    mxAssert(nrows >= blkstart[nblk], "X size mismatch");
    X.jc = mxGetJc(X_IN);
    X.ir = mxGetIr(X_IN);
/* ------------------------------------------------------------
   Get XqjcPr, pointing to start of Lorentz blocks in X.
   ------------------------------------------------------------ */
    mxAssert(nrhs >= NPARIN, "ddot with sparse X requires more input arguments.");
    mxAssert(mxGetM(XBLKJC_IN) == m && mxGetN(XBLKJC_IN) >= 3, "Xjc size mismatch");
    XjcPr = mxGetPr(XBLKJC_IN) + m;      /* Point to Xjc(:,2) */
/* ------------------------------------------------------------
   Allocate working arrays:
   mwIndex xjc1(2*m), xblk(qDim).
   ------------------------------------------------------------ */
    xjc1 = (mwIndex *) mxCalloc(MAX(2*m,1), sizeof(mwIndex) );
    xblk = (mwIndex *) mxCalloc(MAX(qDim,1), sizeof(mwIndex) );
/* ------------------------------------------------------------
   Convert double to mwIndex:
   ------------------------------------------------------------ */
    for(i = 0; i < 2*m; i++)
      xjc1[i] = (mwIndex) XjcPr[i];                /* double to mwIndex */
/* ------------------------------------------------------------
   Let k = xblk(j-blkstart[0]) iff
   blkstart[k] <= j < blkstart[k+1], k=0:nblk-1.
   ------------------------------------------------------------ */
    j = 0;
    for(k = 0; k < nblk; k++){
      i = blkstart[k+1] - blkstart[0];
      while(j < i)
        xblk[j++] = k;
    }
/* ------------------------------------------------------------
   Let maxnnz := sum(xjc1(:,2)-xjc1(:,1)).
   Create sparse output ddotX(nblk,m,maxnnz)
   ------------------------------------------------------------ */
    maxnnz = 0;
    for(i = 0; i < m; i++)
      maxnnz += xjc1[m+i] - xjc1[i];
    maxnnz = MAX(1, maxnnz);
    DDOTX_OUT = mxCreateSparse(nblk,m, maxnnz,mxREAL);
    ddotx.jc = mxGetJc(DDOTX_OUT);
    ddotx.ir = mxGetIr(DDOTX_OUT);
    ddotx.pr = mxGetPr(DDOTX_OUT);
/* ------------------------------------------------------------
   The real job:
   ------------------------------------------------------------ */
    spddotxj(ddotx, d, X.ir, X.pr,xjc1,xjc1+m, xblk,blkstart,nblk,m);
/* ------------------------------------------------------------
   REALLOC (shrink) ddotx to ddotx.jc[m] nonzeros.
   ------------------------------------------------------------ */
    maxnnz = MAX(1,ddotx.jc[m]);
    if((ddotx.ir = (mwIndex *) mxRealloc(ddotx.ir, maxnnz * sizeof(mwIndex))) == NULL)
      mexErrMsgTxt("Memory allocation error");
    if((ddotx.pr = (double *) mxRealloc(ddotx.pr, maxnnz*sizeof(double)))
       == NULL)
      mexErrMsgTxt("Memory allocation error");
    mxSetPr(DDOTX_OUT,ddotx.pr);
    mxSetIr(DDOTX_OUT,ddotx.ir);
    mxSetNzmax(DDOTX_OUT,maxnnz);
/* ------------------------------------------------------------
   Release working arrays (SPARSE PART).
   ------------------------------------------------------------ */
    mxFree(xjc1);
    mxFree(xblk);
  }
/* ------------------------------------------------------------
   Release common working arrays.
   ------------------------------------------------------------ */
  mxFree(blkstart);
}
