/*
 x = symbbwblk(L, b)

 L is structure with m x m sparse Cholesky factor L.L, b is m x k sparse rhs.
 Computes sparsity structure of x = L.L' \ b in linear time.
 Property used: L(i,j) != 0, L(k,j) != 0 ==> L(k,i) != 0 for j < i < k.
 Thus, for each supernodal column L(:,j), we only need the first nonzero
 supernode i > j; the remaining ones have already been taken care of
 while backward-solving i.


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

#define X_OUT  plhs[0]
#define NPAROUT 1

#define L_IN      prhs[0]
#define B_IN      prhs[1]
#define NPARIN 2

/* ************************************************************
   PROCEDURE getSnodeBelow  - Get supernodal partition and first
     below-diag supernode of L.
   INPUT
     ljc, lir - uncompressed nz-structure of L
     xsuper - supernodal partition (length nsuper+1).
     nsuper - number of supernodes
   OUTPUT
     snodebelow - length nsuper, the first nonzero supernode below the
       supernodal diag block. Is nsuper if none.
     snode  - length m = xsuper[nsuper+1]. Maps subnode to its supernode.
   ************************************************************ */
void getSnodeBelow(int *snodebelow, int *snode,
                   const int *ljc, const int *lir, const int *xsuper,
                   const int nsuper)
{
  int j, jsup, ix, nj, jcol;
/* ------------------------------------------------------------
   SNODE: map each column to the supernode containing it
   ------------------------------------------------------------ */
  j = xsuper[0];
  for(jsup = 0; jsup < nsuper; jsup++){
    while(j < xsuper[jsup + 1])
      snode[j++] = jsup;
  }
/* ------------------------------------------------------------
   FIRST SUPERNODAL SUBSCRIPT BELOW DIAGONAL:
    Let snodebelow[jsup] = snode[lir[ljc[jcol] + nj]],
    where xsuper[jsup:jsup+1] = jcol : jcol+nj.
   ------------------------------------------------------------ */
  for(jsup = 0, jcol = 0; jsup < nsuper; jsup++){
    nj = xsuper[jsup + 1] - jcol;       /* nj = length snode jsup */
    ix = ljc[jcol] + nj;                /* points just below diag block */
    if(ix < ljc[jcol + 1])              /* Anything below diag-block ? */
      snodebelow[jsup] = snode[lir[ix]];
    else
      snodebelow[jsup] = nsuper;        /* nsuper = nothing below diag */
    jcol += nj;           /* jcol = xsuper[jsup+1] */
  }
}

/* ************************************************************
   PROCEDURE bwnzsuper - Compute sparsity structure of L'\b, by
       determining nonzero-supernodes, and "1 beyond" subnodes within
       them (each supernode is a dense diag block in L).
   INPUT
     bir, bnnz     - bir(bnnz) lists the row-indices of vector b.
     snode, xsuper - supernodal partition.
         xsuper(nsuper+1): xsuper(j):xsuper(j+1)-1 is jth supernode
         snode(m): j=snode(i) means xsuper(j) <= i < xsuper(j+1).
     snodebelow  - First nonzero supernode below jth diagonal block, if any.
        Otherwise, snodebelow[j] = nsuper.
   UPDATED
     processed - char(nsuper+1) array. On input all-0, on output
       processed[jsup] = 1 iff jsup is in nz structure of L'\b.
   OUTPUT
     snodeto - Length nsuper array. snodeto[jsup]-1 is last relevant
       subnode of each supernode jsup where processed[jsup]=1.
   ************************************************************ */
void bwnzsuper(int *snodeto, char *processed,
                const int *bir, const int bnnz,
                const int *snode, const int *xsuper,
                const int *snodebelow)
{
  int inz,jsup,lastsub;
/* ------------------------------------------------------------
   Browse all supernodes jsup = snode[ bir[ inz ] ], letting snodeto[jsup]
   point beyond its last subnode in bir.
   NOTE: if b = e_i, then only jsup=snode[i] has snodeto, the other
   supernodes will be full.
   ------------------------------------------------------------ */
  if(bnnz <= 0)
    return;
  lastsub = bir[bnnz - 1];
  inz = 0;
  while(inz < bnnz){
    jsup = snode[bir[inz]];                  /* jsup = current snode */
    if(xsuper[jsup + 1] <= lastsub)          /* Let inz point to next snode */
      while(bir[++inz] < xsuper[jsup + 1]);
    else
      inz = bnnz;
    snodeto[jsup] = bir[inz - 1] + 1;        /* last subscript + 1 */
    processed[jsup] = 1;
  }
/* ------------------------------------------------------------
   Symbolic supernodal backward solve:
   1) If !processed[snodebelow[jsup]], then no contribution from
     L(:,jsup), since remaining sparsity structure is part of snodebelow[jsup]
   2) Otherwise, the whole supernode gets nonzero, since L(jsup,jsup) is
     supernodal diag block.
   ------------------------------------------------------------ */
  for(--jsup; jsup >= 0; jsup--)
    if(processed[snodebelow[jsup]]){
      processed[jsup] = 1;
      snodeto[jsup] = xsuper[jsup + 1];      /* Include full supernode */
    }
}

/* ************************************************************
   PROCEDURE symbbwmat - Computes nz-structure of x = L'\b.
   INPUT
     bjc, bir - nz-structure of m x n RHS-matrix b.
     snode, xsuper - supernodal partition.
         xsuper(nsuper+1): xsuper(j):xsuper(j+1)-1 is jth supernode
         snode(m): j=snode(i) means xsuper(j) <= i < xsuper(j+1).
     xlindx,lindx  - compressed subscript array.
         xlindx(nsuper+1): lindx(xlindx(j):xlindx(j+1)-1) are the subscripts
         for supernode j.
     nsuper   - number of supernodes
     m,n      - size(b), m rows, n columns.
   OUTPUT
     xjc     - n+1-array: Start of each column in *pxir.
     *pxir   - length *pmaxnnz array of row-indices.
   UPDATED
     *pmaxnnz - The allocated number of entries in *pxir. Will be changed
        by this function to the exact number needed (but s.t. maxnnz >= 1).
   WORK
     snodeto - int(nsuper)
     processed - char(nsuper+1)
   ************************************************************ */
void symbbwmat(int *xjc, int **pxir,int *pmaxnnz,
               const int *bjc, const int *bir,
               const int *snode, const int *xsuper,
               const int *snodebelow,
               const int nsuper, const int m, const int n,
               int *snodeto, char *processed)
{
  int i,j,k,inz, maxnnz;
  int *xir;
/* ------------------------------------------------------------
   INIT: processed = 0, xir = *pxir, maxnnz = *pmaxnnz.
   processes[nsuper]=0 is used for "empty"-flag "snodebelow[j]=nsuper".
   ------------------------------------------------------------ */
  memset(processed, 0, (nsuper + 1) * sizeof(char));
  xir = *pxir;
  maxnnz = *pmaxnnz;
/* ------------------------------------------------------------
   For each column j, compute nz-structure of L'\b(:,j) into xir.
   First make sure that xir has enough (at least m) available entries.
   ------------------------------------------------------------ */
  inz = 0;
  for(j = 0; j < n; j++){
    xjc[j] = inz;
    if(inz + m > maxnnz){
      maxnnz += inz + m;        /* required + old amount */
      xir = (int *) mxRealloc(xir, maxnnz * sizeof(int));
    }
/* ------------------------------------------------------------
   Find all nz-supernodes in L'\b(:,j).
   ------------------------------------------------------------ */
    bwnzsuper(snodeto, processed, bir + bjc[j], bjc[j+1]-bjc[j],
              snode, xsuper, snodebelow);
/* ------------------------------------------------------------
   For each nz-supernode, write the row-indices "< snodeto" into xir.
   ------------------------------------------------------------ */
    for(k = 0; k < nsuper; k++)
      if(processed[k]){
        processed[k] = 0;
        for(i = xsuper[k]; i < snodeto[k]; i++)
          xir[inz++] = i;
      }
  }
/* ------------------------------------------------------------
   FINALLY: close last column in xir, Realloc xir to the actual
   maxnnz:= xjc[n], and return.
   ------------------------------------------------------------ */
  xjc[n] = inz;
  if(inz < maxnnz){
    maxnnz = MAX(inz,1);           /* avoid realloc to NULL */
    xir = (int *) mxRealloc(xir, maxnnz * sizeof(int));
  }
  *pxir = xir;
  *pmaxnnz = maxnnz;
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  const mxArray *L_FIELD;
  int maxnnz, i,j, nsuper,m,n;
  const int *ljc,*lir,*bjc,*bir;
  int *xjc,*xir, *snode,*snodebelow, *iwork,*xsuper;
  char *cwork;
  double *xpr;
  const double *xsuperPr;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "symbfwblk requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "symbfwblk produces 1 output argument");
/* ------------------------------------------------------------
   Get rhs-input B
   ------------------------------------------------------------ */
  mxAssert(mxIsSparse(B_IN), "B must be sparse");
  m = mxGetM(B_IN);
  n = mxGetN(B_IN);
  bjc = mxGetJc(B_IN);
  bir = mxGetIr(B_IN);
/* ------------------------------------------------------------
   Disassemble block Cholesky structure L
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(L_IN), "Parameter `L' should be a structure.");
  L_FIELD = mxGetField(L_IN,0,"L"); 
  mxAssert( L_FIELD != NULL, "Missing field L.L.");           /* L.L */
  mxAssert( m == mxGetM(L_FIELD) && m == mxGetN(L_FIELD), "Size L.L mismatch.");
  mxAssert(mxIsSparse(L_FIELD), "L.L should be sparse.");
  ljc = mxGetJc(L_FIELD);
  lir = mxGetIr(L_FIELD);
  L_FIELD = mxGetField(L_IN,0,"xsuper"); 
  mxAssert( L_FIELD != NULL, "Missing field L.xsuper.");      /* L.xsuper */
  nsuper = mxGetM(L_FIELD) * mxGetN(L_FIELD) - 1;
  mxAssert( nsuper <= m, "Size L.xsuper mismatch.");
  xsuperPr = mxGetPr(L_FIELD);
/* ------------------------------------------------------------
   Allocate int-part of sparse output matrix X(m x n)
   Heuristically set nnz to nnz(B) + 4*m.
   ------------------------------------------------------------ */
  maxnnz = bjc[n] + 4 * m;
  xjc = (int *) mxCalloc(n + 1, sizeof(int));
  xir = (int *) mxCalloc(maxnnz, sizeof(int));
/* ------------------------------------------------------------
   Allocate working arrays:
   int snode(m), xsuper(nsuper+1), snodebelow(nsuper),
   iwork(nsuper).
   char cwork(nsuper+1).
   ------------------------------------------------------------ */
  snode     = (int *) mxCalloc(m,sizeof(int)); 
  xsuper    = (int *) mxCalloc(nsuper+1,sizeof(int));
  snodebelow = (int *) mxCalloc(nsuper,sizeof(int));
  iwork = (int *) mxCalloc(nsuper, sizeof(int));
  cwork = (char *) mxCalloc(nsuper+1, sizeof(char));
/* ------------------------------------------------------------
   Convert XSUPER to integer and C-Style
   ------------------------------------------------------------ */
  for(i = 0; i <= nsuper; i++){
    j =  xsuperPr[i];
    xsuper[i] = --j;
  }
/* ------------------------------------------------------------
   Create "snode" from xsuper, and get "first-below-diag" 
   supernodal subscript snodebelow (snodebelow[j]==nsuper means none).
   This is enough to determine the nz-pattern of the backward-solve.
   ------------------------------------------------------------ */
  getSnodeBelow(snodebelow,snode, ljc,lir,xsuper,nsuper);
/* ------------------------------------------------------------
   Compute nz structure after backward solve
   ------------------------------------------------------------ */
  symbbwmat(xjc, &xir, &maxnnz, bjc, bir, snode, xsuper,
            snodebelow, nsuper, m, n, iwork, cwork);
/* ------------------------------------------------------------
   Create output matrix x
   ------------------------------------------------------------ */
  X_OUT = mxCreateSparse(m,n, 1,mxREAL);
  mxFree(mxGetJc(X_OUT));                    /* jc */
  mxFree(mxGetIr(X_OUT));                    /* ir */
  mxFree(mxGetPr(X_OUT));                    /* pr */
  xpr = (double *) mxCalloc(maxnnz,sizeof(double));
  mxSetJc(X_OUT, xjc);
  mxSetIr(X_OUT, xir);
  mxSetPr(X_OUT, xpr);
  mxSetNzmax(X_OUT, maxnnz);
  for(i = 0; i < maxnnz; i++)
    xpr[i] = 1.0;
/* ------------------------------------------------------------
   Release working arrays.
   ------------------------------------------------------------ */
  mxFree(cwork);
  mxFree(iwork);
  mxFree(snodebelow);
  mxFree(xsuper);
  mxFree(snode);
}
