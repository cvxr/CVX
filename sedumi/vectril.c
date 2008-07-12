/*
   y = vectril(x,K)
   For the PSD submatrices, we let Yk = tril(Xk,0) + triu(Xk,1)'
   Complex numbers are stored in SeDuMi's (real Xk; imag Xk) format.

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

#include <math.h>
#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define Y_OUT plhs[0]
#define NPAROUT 1

#define X_IN prhs[0]
#define K_IN prhs[1]
#define NPARIN 2

/* ************************************************************
   PROCEDURE sptriujcT - Compute starting positions for storing
     nonzeros of rows triu(X(i,:),1). Note: X is nxn block stored
     as subvector in a long vector.
   INPUT
     xir   - integer subscript array of length xjc1.
     xjc0,xjc1 - xjc0 points to n x n sparse matrix X in xir. 
     first - subscript of X(0,0); vec(X) is stored as subvector in xir.
     n - order of X.
   OUTPUT
     triujc - length n-1 integer array; triujc[i] is starting index
       for storing triu(X(i+1,:),1), i=0:n-2.
       It starts at triujc[0] = nnz(tril(X)).
   RETURNS nnz(X). Note that
   nnz(X) = triujc[n-2] + nnz(triu(X(n-1,:),1)) <= xjc1 - xjc0.
   ************************************************************ */
mwIndex sptriujcT(mwIndex *triujc, const mwIndex *xir, const mwIndex xjc0, const mwIndex xjc1,
              const mwIndex first,const mwIndex n)
{
  mwIndex i,j,inz,jnz,jfirst,jlast;
/* ------------------------------------------------------------
   Observe that triu(X(n,:))=[] and hence only use triujc[0:n-2].
   ------------------------------------------------------------ */
  mxAssert(n > 0,"");
  memset(triujc, 0, (n-1) * sizeof(mwIndex));
  jnz = 0;        /* jnz = nnz(tril(X)) */
  jlast = 0;      /* index right after last activated column */
  for(inz = xjc0; inz < xjc1; inz++){
/* ------------------------------------------------------------
   Translate vec(x)-subscript into (i,j) subscript
   ------------------------------------------------------------ */
    if((i = xir[inz]) >= jlast){      /* move to new z-column */
      j = (i-first) / n;              /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all columns of X */
      jfirst = first + n * j;
      jlast = jfirst + n;
    }
    mxAssert(i >= jfirst,"");
    i -= jfirst;
/* ------------------------------------------------------------
   x_{i,j} with i < j --> count as nonzero in row i of triu(x)
   ------------------------------------------------------------ */
    if(i < j)              /* Upper triangular nonzero in row i */
      ++triujc[i];    /* 0 <= i <= j-1 <= n-2 */
/* ------------------------------------------------------------
   x_{i,j} with i >= j -->  count as nonzero in tril(x).
   ------------------------------------------------------------ */
    else{
      ++jnz;
    }
  }
/* ------------------------------------------------------------
   Now jnz = nnz(tril(X)) and triujc[i] = nnz(triu(X(i+1,:),1)).
   We will now let triujc point to 1st nonzero position for storing
   jth row of triu(X,1) in row-wise format, as follows:
   [jnz, jnz + triujc[0], .., jnz + triujc[n-2]].
   ------------------------------------------------------------ */
  for(i = 0; i < n-1; i++){
    j = triujc[i];        /* nnz(triu(X(i,:),1)) */
    triujc[i] = jnz;      /* position to store triu(X(i,:),1) */
    jnz += j;            /* point to next available position */
  }
/* ------------------------------------------------------------
   Now jnz = triujc[n-1, NEW] + triujc[n-1,OLD] = nnz(X).
   ------------------------------------------------------------ */
  mxAssert(jnz == inz - xjc0,"");
  return jnz;
}

/* ************************************************************
   PROCEDURE sptrilandtriu - Compute y = [tril(X); triu(X,1)']
   INPUT
     xir   - integer subscript array of length xjc1.
     xpr   - vector of xjc1 nonzeros.
     xjc0,xjc1 - xjc0 points to n x n sparse matrix X in xir. 
     first - subscript of X(0,0); vec(X) is stored as subvector in xir.
     n - order of X.
   UPDATED
     triujc - length n-1 integer array. INPUT: triujc[i] is starting index
       for storing triu(X(i+1,:),1), i=0:n-2. It starts at triujc[0] =
       nnz(tril(X)). OUTPUT: triujc[i] points beyond last nonzero
       written in triu(X(i,:),1)'. Thus, triujc[n-2] = nnz(y).
   OUTPUT
     yir, ypr - length nnz(y) sparse vector, y = [tril(X); triu(X,1)'].
   RETURNS tilxnnz = nnz(tril(X)) if !skew, and nnz(tril(X,-1)) if skew.
   ************************************************************ */
mwIndex sptrilandtriu(mwIndex *yir, double *ypr, mwIndex *triujc, const mwIndex *xir,
                   const double *xpr, const mwIndex xjc0, const mwIndex xjc1,
                   const mwIndex first, const mwIndex n, bool skew)
{
  mwIndex i,j,inz,jnz,knz,jfirst,jlast;
/* ------------------------------------------------------------
   Store [tril(X); triu(X,1)'] into y,
   using the column pointers of triu(X,1)' as computed in triujc.
   ------------------------------------------------------------ */
  jnz = 0;
  jlast = 0;      /* index right after last activated column */
  for(inz = xjc0; inz < xjc1; inz++){
/* ------------------------------------------------------------
   Translate vec(x)-subscript into (i,j) subscript
   ------------------------------------------------------------ */
    if((i = xir[inz]) >= jlast){      /* move to new z-column */
      j = (i-first) / n;              /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all columns of X */
      jfirst = first + n * j;
      jlast = jfirst + n;
    }
    i -= jfirst;
/* ------------------------------------------------------------
   x_{i,j} with i < j --> store nonzero in row i of triu(x)
   ------------------------------------------------------------ */
    if(i < j){              /* Upper triangular nonzero in row i */
      knz = triujc[i];
      ++triujc[i];
      yir[knz] = first + n * i + j;   /* vec-index  after transposing */
      ypr[knz] = xpr[inz];
    }
/* ------------------------------------------------------------
   x_{i,j} with i >= j --> store nonzero of tril(x) or, for skew, tril(x,-1).
   ------------------------------------------------------------ */
    else if(!skew || i > j){
      yir[jnz] = jfirst + i;   /* vec-index */
      ypr[jnz] = xpr[inz];
      ++jnz;
    }
  }
/* ------------------------------------------------------------
   Return number of tril-nonzeros written.
   !skew ==> nnz(tril(X));  skew ==> nnz(tril(X,-1)).
   ------------------------------------------------------------ */
  return jnz;
}

/* ************************************************************
   PROCEDURE spadd - Let z = x + y
   INPUT
     xir, xpr, xnnz - sparse vector
     yir, ypr, ynnz - sparse vector
     iwsize - ynnz+2 + floor(log_2(ynnz+1))
   OUTPUT
     zir, zpr - length znnz arrays containing sparse output z = x + y.
   WORK
     iwork - length iwsize working array
   RETURNS znnz
   ************************************************************ */
mwIndex spadd2(mwIndex *zir, double *zpr, const mwIndex *xir, const double *xpr,
           const mwIndex xnnz, const mwIndex *yir, const double *ypr, const mwIndex ynnz,
           mwIndex iwsize, bool *cfound, mwIndex *iwork)
{
  mwIndex inz,jnz,knz, i;
  mwIndex * yinx;
/* ------------------------------------------------------------
   Partition working array [yinx(ynnz+2), iwork(log_2(ynnz+1))].
   ------------------------------------------------------------ */
  yinx = iwork;
  iwork += ynnz + 2;
  iwsize -= ynnz + 2;
  intmbsearch(yinx, cfound, xir, xnnz, yir, ynnz, iwork, iwsize);
  jnz = yinx[1];
  memcpy(zir, xir, jnz * sizeof(mwIndex));
  memcpy(zpr, xpr, jnz * sizeof(double));
  for(i = 0; i < ynnz; i++){
    inz = yinx[i+1];
    if(cfound[i])
      zpr[jnz] = xpr[inz++] + ypr[i];
    else
      zpr[jnz] = ypr[i];
    zir[jnz++] = yir[i];
    knz = yinx[i+2]-inz;
    memcpy(zir + jnz, xir + inz, knz * sizeof(mwIndex));
    memcpy(zpr + jnz, xpr + inz, knz * sizeof(double));
    jnz += knz;
  }
  return jnz;
}
      
mwIndex spadd(mwIndex *zir, double *zpr, const mwIndex *xir, const double *xpr,
          const mwIndex xnnz, const mwIndex *yir, const double *ypr, const mwIndex ynnz,
          const mwIndex iwsize, bool *cfound, mwIndex *iwork)
{
  if(xnnz < ynnz)
    return spadd2(zir,zpr, yir,ypr,ynnz, xir,xpr,xnnz, iwsize, cfound, iwork);
  else
    return spadd2(zir,zpr, xir,xpr,xnnz, yir,ypr,ynnz, iwsize, cfound, iwork);
}

/* ************************************************************
   PROCEDURE spsub - Let z = x - y
   INPUT
     xir, xpr, xnnz - sparse vector
     yir, ypr, ynnz - sparse vector
     iwsize - ynnz+2 + floor(log_2(ynnz+1))
   OUTPUT
     zir, zpr - length znnz arrays containing sparse output z = x - y.
   WORK
     iwork - length iwsize working array
   RETURNS znnz
   ************************************************************ */
mwIndex spsub(mwIndex *zir, double *zpr, const mwIndex *xir, const double *xpr,
          const mwIndex xnnz, const mwIndex *yir, const double *ypr, const mwIndex ynnz,
          mwIndex iwsize, bool *cfound, mwIndex *iwork)
{
  mwIndex inz,jnz,knz, i;
  mwIndex *yinx;
/* ------------------------------------------------------------
   Partition working array [yinx(ynnz+2), iwork(log_2(ynnz+1))].
   ------------------------------------------------------------ */
  yinx = iwork;
  iwork += ynnz + 2;
  iwsize -= ynnz + 2;
  intmbsearch(yinx, cfound, xir, xnnz, yir, ynnz, iwork, iwsize);
  jnz = yinx[1];
  memcpy(zir, xir, jnz * sizeof(mwIndex));
  memcpy(zpr, xpr, jnz * sizeof(double));
  for(i = 0; i < ynnz; i++){
    inz = yinx[i+1];
    if(cfound[i])
      zpr[jnz] = xpr[inz++] - ypr[i];
    else
      zpr[jnz] = -ypr[i];
    zir[jnz++] = yir[i];
    knz = yinx[i+2]-inz;
    memcpy(zir + jnz, xir + inz, knz * sizeof(mwIndex));
    memcpy(zpr + jnz, xpr + inz, knz * sizeof(double));
    jnz += knz;
  }
  return jnz;
}
      
/* ************************************************************
   PROCEDURE sptotril - For sparse x=vec(X), lets
     z = vec( tril(X) + triu(X,1)' ). If skew = 1 then
     z = vec( tril(X) - triu(X)' ).
   INPUT
     xir, xpr, pxjc0, xjc1 - sparse input vector, *pxjc0 points to
       first nonzero of vectorized matrix X.
     first - subscript of X(1,1) in long vector x.
     n - order of n x n matrix X.
     skew - if 1, then set SUBTRACT triu(X,1)' and set diag(z)=all-0.
     iwsize - n + xnnz + 1+nnz(triu(X,1)) + log_2(1+nnz(triu(X,1))).
        Observe that nnz(triu(X,1)) <= MIN(n*(n-1)/2, xnnz), and
        xnnz <= MIN(n^2, xjc1-*pxjc0). Thus
        iwsize <= n*(2*n+1)+log_2(1+n*(n-1)/2).
   OUTPUT
     zir - length znnz mwIndex array, subscripts of z := vec(tril(x)+triu(x,1)').
     zpr - length znnz vector, nonzeros of z.
   WORK
     cwork - length nnz(triu(X,1)) <= n*(n-1)/2 char array.
     iwork - length iwsize integer working array
     ypr - length xnnz vector; xnnz <= n^2.
   RETURNS znnz
   ************************************************************ */
mwIndex sptotril(mwIndex *zir, double *zpr, const mwIndex *xir, const double *xpr,
             mwIndex *pxjc0, const mwIndex xjc1, const mwIndex first, const mwIndex n,
             bool skew, mwIndex iwsize, bool *cwork, mwIndex *iwork,
             double *ypr)
{
  mwIndex xjc0, xnnz, trilnnz, triujc0;
  mwIndex *triujc, *yir;
/* ------------------------------------------------------------
   Let iwork[0:n-2] point to row-starts for storing triu(X,1)
   row-wise. Let xnnz be nnz(X). Update *pxjc0 to point beyond this
   block
   ------------------------------------------------------------ */
  xjc0 = *pxjc0;
  xnnz = sptriujcT(iwork, xir, xjc0, xjc1, first, n);
  *pxjc0 = xjc0 + xnnz;
/* ------------------------------------------------------------
   Partition integer working array
   ------------------------------------------------------------ */
  triujc = iwork;
  yir = iwork + (n-1);
  iwork = yir + xnnz;
  iwsize -= n-1 + xnnz;
/* ------------------------------------------------------------
   ------------------------------------------------------------ */
  if(n > 1)
    triujc0 = triujc[0];
  else
    triujc0 = xnnz;            /* 1 x 1 matrix --> triu(X)=[] */
  trilnnz = sptrilandtriu(yir, ypr, triujc, xir,xpr,xjc0,xjc1, first,n, skew);
  if(!skew)
    return spadd(zir,zpr, yir,ypr, trilnnz, yir+triujc0,ypr+triujc0,
                 xnnz-triujc0, iwsize, cwork, iwork);
  else
    return spsub(zir,zpr, yir,ypr, trilnnz, yir+triujc0,ypr+triujc0,
                 xnnz-triujc0, iwsize, cwork, iwork);
}

/* ************************************************************
   PROCEDURE vectril - Applies sptotril(xk) for each PSD block k.
       On output, each PSD block is lower triangular, i.e.
       Zk = tril(Xk+Xk')/2.
   INPUT
     xir,xpr,xnnz
     psdNL - K.s
     blkstart - length psdN+1 array. PSD block k has subscripts
       blkstart[k]:blkstart[k+1]-1.
     isblk - length psdDim array, with k = xblk(i-blkstart[0]) iff
       blkstart[k] <= i < blkstart[k+1], k=0:psdN-1.
     rpsdN - number of real PSD blocks
     psdN - number of PSD blocks
     iwsize - maxn*(2*maxn+1)+log_2(1+maxn*(maxn-1)/2), where maxn := max(K.s).
   OUTPUT
     zir - length znnz <= xnnz mwIndex array: subscripts of z = vectril(x).
     zpr - length znnz <= xnnz vector: nonzeros of z = vectril(x).
   WORKING ARRAYS
     cwork - length maxn*(maxn-1)/2 char array, where maxn := max(K.s).
     iwork - length iwsize integer working array
     fwork - length max(K.s.^2) vector. (Note: not double for Hermitian
       blocks, since we treat real and imag parts seperately.)
   RETURNS znnz
   ************************************************************ */
mwIndex vectril(mwIndex *zir, double *zpr, const mwIndex *xir, const double *xpr,
            const mwIndex xnnz, const mwIndex *psdNL,
            const mwIndex *blkstart, const mwIndex *isblk,
            const mwIndex rpsdN, const mwIndex psdN, const mwIndex iwsize,
            bool *cwork, mwIndex *iwork, double *fwork)
{
  mwIndex inz, jnz, k, nk;
/* ------------------------------------------------------------
   Copy f,l,q,r parts without change. Let inz point to first
   PSD-nonzero in x, jnz in z.
   ------------------------------------------------------------ */
  inz = 0;                                   /* pointer into x */
  intbsearch(&inz, xir, xnnz, blkstart[0]);  /* inz points to start PSD */
  isblk -= blkstart[0];
  memcpy(zir, xir, inz * sizeof(mwIndex));
  memcpy(zpr, xpr, inz * sizeof(double));
  jnz = inz;                        /* jnz points to start PSD in z */
/* ------------------------------------------------------------
   Process all PSD blocks
   ------------------------------------------------------------ */
  while(inz < xnnz){
    k = isblk[xir[inz]];
    nk = psdNL[k];
    jnz += sptotril(zir + jnz, zpr + jnz, xir, xpr, &inz, xnnz, blkstart[k],
                    nk,(bool)0, iwsize, cwork, iwork, fwork);
/* ------------------------------------------------------------
   For the imaginary part, we do a skew transpose: tril(IM Xk)-triu(IM Xk)'.
   This will make the diagonal of the imaginary block zero.
   ------------------------------------------------------------ */
    if(k >= rpsdN){
      jnz += sptotril(zir + jnz, zpr + jnz, xir, xpr, &inz, xnnz,
                      blkstart[k]+SQR(nk), nk,(bool)1, iwsize, cwork, iwork, fwork);
    }
  }
  return jnz;
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
     y = vectril(x,K)

   For the PSD submatrices, we let Yk = tril(Xk+Xk').
   Complex numbers are stored as vec([real(Xk) imag(Xk)]).
   NB: x and y are sparse.
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  mwIndex i, j, k, jnz, m,lenfull, firstPSD, maxn, iwsize;
  jcir x,y;
  mwIndex *iwork, *psdNL, *blkstart, *xblk;
  bool *cwork;
  double *fwork;
  coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "vectril requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "vectril produces less output arguments");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Compute statistics based on cone K structure
   ------------------------------------------------------------ */
  firstPSD = cK.frN + cK.lpN + cK.qDim;
  for(i = 0; i < cK.rconeN; i++)        /* add dim of rotated cone */
    firstPSD += (mwIndex) cK.rconeNL[i];
  lenfull =  firstPSD + cK.rDim + cK.hDim;
/* ------------------------------------------------------------
   Get inputs x, blkstart
   ------------------------------------------------------------ */
  mxAssert(mxGetM(X_IN) == lenfull, "X size mismatch.");
  m = mxGetN(X_IN);                       /* number of columns to handle */
  mxAssert( mxIsSparse(X_IN), "X should be sparse.");
  x.pr = mxGetPr(X_IN);
  x.jc = mxGetJc(X_IN);
  x.ir = mxGetIr(X_IN);
/* ------------------------------------------------------------
   Allocate output Y = sparse([],[],[],length(x),m,nnz(x))
   ------------------------------------------------------------ */
  Y_OUT = mxCreateSparse(lenfull, m, x.jc[m], mxREAL);
  y.pr = mxGetPr(Y_OUT);
  y.jc = mxGetJc(Y_OUT);
  y.ir = mxGetIr(Y_OUT);
  y.jc[0] = 0;
/* ------------------------------------------------------------
   If x = [], then we are ready with y=[]. Otherwise, proceed:
   ------------------------------------------------------------ */
  if(x.jc[m] > 0){
/* ------------------------------------------------------------
   Allocate iwork[iwsize],
   iwsize := maxn*(2*maxn+1)+log_2(1+maxn*(maxn-1)/2), where maxn := max(K.s);
   cwork[maxn*(maxn-1)/2], fwork(maxn^2), mwIndex psdNL(length(K.s)).
   mwIndex blkstart(sdpN+1), xblk(sdpDim)
   ------------------------------------------------------------ */
    maxn = MAX(cK.rMaxn,cK.hMaxn);
    iwsize = (mwIndex) (log(1.0 + maxn*(maxn-1)/2) / log(2.0));
    iwsize += maxn * (2*maxn+1);
    iwork = (mwIndex *) mxCalloc(MAX(1,iwsize), sizeof(mwIndex));
    cwork = (bool *) mxCalloc(MAX(1,maxn*(maxn-1)/2), sizeof(bool));
    fwork = (double *) mxCalloc(MAX(1,SQR(maxn)), sizeof(double));
    psdNL = (mwIndex *) mxCalloc(MAX(1,cK.sdpN), sizeof(mwIndex));
    blkstart = (mwIndex *) mxCalloc(1 + cK.sdpN, sizeof(mwIndex));
    xblk = (mwIndex *) mxCalloc(MAX(1,cK.rDim + cK.hDim), sizeof(mwIndex));
/* ------------------------------------------------------------
   double -> mwIndex for K.s
   ------------------------------------------------------------ */
    for(i = 0; i < cK.sdpN; i++)
      psdNL[i] = (mwIndex) cK.sdpNL[i];
/* ------------------------------------------------------------
   Let k = xblk(j-blkstart[0]) iff
   blkstart[k] <= j < blkstart[k+1], k=0:psdN-1.
   ------------------------------------------------------------ */
    j = firstPSD;
    for(i = 0; i < cK.rsdpN; i++){     /* real sym */
      blkstart[i] = j;
      j += SQR(psdNL[i]);
    }
    for(; i < cK.sdpN; i++){            /* complex herm. */
      blkstart[i] = j;
      j += 2*SQR(psdNL[i]);
    }
    blkstart[cK.sdpN] = j;
    mxAssert(j - firstPSD == cK.rDim + cK.hDim, "Size mismatch blkstart, K.");
    j = 0;
    for(k = 0; k < cK.sdpN; k++){
      i = blkstart[k+1] - blkstart[0];
      while(j < i)
        xblk[j++] = k;
    }
/* ------------------------------------------------------------
   Let y(:,i)= vectril(x(:,i)), for i=1:m.
   ------------------------------------------------------------ */
    jnz = 0;            /* points into y */
    for(i = 0; i < m; i++){
      y.jc[i] = jnz;
      jnz += vectril(y.ir+jnz,y.pr+jnz, x.ir+x.jc[i],x.pr+x.jc[i],
                     x.jc[i+1]-x.jc[i],
                     psdNL, blkstart, xblk, cK.rsdpN,cK.sdpN, iwsize,
                     cwork, iwork, fwork);
    }
    y.jc[m] = jnz;    /* nnz written into y */
    mxAssert(jnz <= x.jc[m],"");
/* ------------------------------------------------------------
   REALLOC: Shrink Y to its current size
   ------------------------------------------------------------ */
    jnz = MAX(jnz,1);
    if( (y.pr = (double *) mxRealloc(y.pr, jnz*sizeof(double))) == NULL)
      mexErrMsgTxt("Memory reallocation error");
    mxSetPr(Y_OUT,y.pr);
    if( (y.ir = (mwIndex *) mxRealloc(y.ir, jnz*sizeof(mwIndex))) == NULL)
      mexErrMsgTxt("Memory reallocation error");
    mxSetIr(Y_OUT,y.ir);
    mxSetNzmax(Y_OUT,jnz);
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
    mxFree(xblk);
    mxFree(blkstart);
    mxFree(psdNL);
    mxFree(fwork);
    mxFree(iwork);
    mxFree(cwork);
  }
}
