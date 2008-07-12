/*
 L.split = cholsplit(L, cachesiz)

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

#define SPLIT_OUT plhs[0]
#define NPAROUT 1

#define L_IN prhs[0]       /* symbolic Cholesky structure: {L.L, L.xsuper} */
#define CACHESIZ_IN prhs[1]       /* number of KBs cache, e.g. 256 */
#define NPARIN 2

#include "mex.h"
#include <math.h>

/* ************************************************************
   PROCEDURE getsplit - Compute splitting of supernodes, such that
     the columns of the splitted supernode fit entirely into the
     computer-cache (whenever possible). Will reduce memory-retrieval
     time in BLKCHOL2.
   INPUT
     ljc,lir - sparsity structure of m x m matrix L (not compressed).
     xsuper,nsuper - supernodal partition of nodes 1:m.
     cachesiz - should be approx. 90% of cache-size in *DOUBLES* (not KBs).
   OUTPUT
     split - length m integer array. For the start of each splitted
       supernode j. Thus, 1 <= split[j] < xsuper[snode[j] + 1].
       For intermediate nodes, split[j] = 0.
   ************************************************************ */
void getsplit(mwIndex *split, const mwIndex *ljc,const mwIndex *lir,const mwIndex *xsuper,
              const mwIndex nsuper, const mwIndex cachesiz)
{
  mwIndex j,k,ksup,mk, used,nextk;

/* ------------------------------------------------------------
   For each supernode ksup = 1:nsuper, column k=1:m.
   ------------------------------------------------------------ */
  k = 0;
  for(ksup = 0; ksup < nsuper; ksup++){
    mk = ljc[k+1] - ljc[k];                /* length of column k */
    used = 2 * mk;                         /* 1st col counts twice */
    nextk = xsuper[ksup + 1];
    j = k;                                 /* j = start cache-group */
/* ------------------------------------------------------------
   If 1st column in ksup is too big to fit into cache, then pick
   columns together until we arrive at a column k that fits.
   ------------------------------------------------------------ */
    if(used > cachesiz){
      k = j + (used - cachesiz) / 2;       /* k is 1st col that fits in cache*/
      if(k >= nextk)
        k = nextk;                         /* all cols in ksup too long */
      else{
        mk -= k-j;
        used = 2*mk;                       /* start new cache-group at k */
      }
      split[j] = k-j;
      j = k;
    }
    else{
      j = k;
      k++;
      mk--;
    }
/* ------------------------------------------------------------
   Split remainder of supernode into cache-groups that fit into cache.
   ------------------------------------------------------------ */
    for(;k < nextk; k++, mk--)
      if(used + mk < cachesiz)
        used += mk;                    /* add into current cache-group */
      else{
        split[j] = k-j;                /* close cache group at previous col. */
        j = k;                         /* start new cache group */
        used = 2 * mk;                 /* insert its first column */
      }
/* ------------------------------------------------------------
   Close last cache-group in this supernode.
   ------------------------------------------------------------ */
    if(j < nextk){
      split[j] = nextk-j;              /* last cache group in this supernode */
    }
  }
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
  const mxArray *L_FIELD;
  mwIndex i,j, nsuper,m, cachesiz;
  const mwIndex *ljc,*lir;
  mwIndex *xsuper, *split;
  const double *xsuperPr;
  double *splitPr;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "cholsplit requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "cholsplit produces less output arguments.");
/* ------------------------------------------------------------
   Get cachesiz, and transform from KBs into 90% of FLOATS.
   ------------------------------------------------------------ */
  cachesiz = (mwIndex) floor(0.9 * (1024 / sizeof(double)) * mxGetScalar(CACHESIZ_IN));
/* ------------------------------------------------------------
   Disassemble block Cholesky structure L
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(L_IN), "Parameter `L' should be a structure.");
  L_FIELD = mxGetField(L_IN,(mwIndex)0,"L");        /* L.L */
  mxAssert( L_FIELD != NULL, "Missing field L.L.");
  m = mxGetM(L_FIELD);
  mxAssert(m == mxGetN(L_FIELD), "L.L must be square.");
  mxAssert(mxIsSparse(L_FIELD), "L.L should be sparse.");
  ljc = mxGetJc(L_FIELD);
  lir = mxGetIr(L_FIELD);
  L_FIELD = mxGetField(L_IN,(mwIndex)0,"xsuper");          /* L.xsuper */
  mxAssert( L_FIELD != NULL, "Missing field L.xsuper.");
  nsuper = mxGetM(L_FIELD) * mxGetN(L_FIELD) - 1;
  mxAssert( nsuper <= m, "Size L.xsuper mismatch.");
  xsuperPr = mxGetPr(L_FIELD);
/* ------------------------------------------------------------
   Allocate working arrays:
   ------------------------------------------------------------ */
  xsuper    = (mwIndex *) mxCalloc(nsuper+1,sizeof(mwIndex));
  split     = (mwIndex *) mxCalloc(m,sizeof(mwIndex));
/* ------------------------------------------------------------
   Convert XSUPER to integer and C-Style
   ------------------------------------------------------------ */
  for(i = 0; i <= nsuper; i++){
    j =  (mwIndex) xsuperPr[i];
    mxAssert(j>0,"");
    xsuper[i] = --j;
  }
/* ------------------------------------------------------------
   The main job: compute (upper bound on) blkchol-split.
   ------------------------------------------------------------ */
  getsplit(split, ljc,lir,xsuper,nsuper, cachesiz);
/* ------------------------------------------------------------
   create OUTPUT variable SPLIT(m)
   ------------------------------------------------------------ */
  SPLIT_OUT = mxCreateDoubleMatrix(m, (mwSize)1, mxREAL);          /* L.split */
  splitPr = mxGetPr(SPLIT_OUT);
  for(i = 0; i < m; i += j){
    j = split[i];
    splitPr[i] = (double) j;
  }
/* ------------------------------------------------------------
   Release working arrays.
   ------------------------------------------------------------ */
  mxFree(split);
  mxFree(xsuper);
}
