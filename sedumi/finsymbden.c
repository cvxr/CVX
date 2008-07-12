/*
%                                   Lden = finsymbden(LAD,perm,dz,firstq)
% FINSYMBDEN  Updates perm and dz by inserting the
%  last Lorentz trace columns (last columns of LAD). It creates the fields
%  Lden.sign  - +1 for "normal" columns, -1 for Lorentz trace columns
%  Lden.first - First pivot column that will affect this one
%  NOTE: sign and first correspond to columns in LAD (without perm-reordering).
%
% SEE ALSO incorder
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
function Lden = finsymbden(LAD,perm,dz,firstq)

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

#define LDEN_OUT  plhs[0]
#define NPAROUT 1

#define LAD_IN    prhs[0]
#define PERM_IN   prhs[1]
#define DZ_IN     prhs[2]
#define FIRSTQ_IN prhs[3]
#define NPARIN 4

/* ************************************************************
   PROCEDURE getfirstpiv - Find first affecting pivot on column
     j = 1:n. These column numbers are in NON-PIVOTED ORDER, i.e.
     the order in which they appear in Xjc.
   INPUT
     invperm - length n array, yields position in list of nonzeros "dzir".
     xsuper - length n (though it may have n+1) array, partitioning
       of permuted subscripts, is "dzjc".
     Xjc - length n+1 array
     Xir - length Xjc[n] array
     n - number of columns in X.
   OUTPUT
     firstpiv - length n array
   ************************************************************ */
void getfirstpiv(mwIndex *firstpiv, const mwIndex *invperm, const mwIndex *xsuper,
                 const mwIndex *Xjc, const mwIndex *Xir, const mwIndex n)
{
  mwIndex i,j,inz,firstj;
  inz = Xjc[0];                 /* typically inz = 0*/
  for(j = 0; j < n; j++){
/* ------------------------------------------------------------
   Let firstj = min(invperm(find(X(:,j))))
   ------------------------------------------------------------ */
    if(inz < Xjc[j+1]){
      firstj = invperm[Xir[inz]];
      for(++inz; inz < Xjc[j+1]; inz++)
        if((i = invperm[Xir[inz]]) < firstj)
          firstj = i;
/* ------------------------------------------------------------
   First node covering firstj, i.e. xsuper[y] < firstj+1 <= xsuper[y+1],
   with y denoting firstpiv[j].
   ------------------------------------------------------------ */
      firstpiv[j] = 0;             /* search from start */
      intbsearch(firstpiv+j,xsuper+1,n-1,firstj+1);
    }
    else
      firstpiv[j] = n;             /* if all-0 then no affecting pivot */
  }
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
#define NLDEN_FIELDS 4
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  mxArray *LDEN_FIELD;
  mwIndex i,inz, j,m,n,nperm, firstQ, lastQ, nnzdz;
  const mwIndex *LADjc, *LADir, *dzJc, *dzIr;
  mwIndex *invdz,*firstpiv,*perm, *dznewJc;
  double *permPr, *firstPr;
  const char *LdenFieldnames[] = {"LAD","perm","dz","first"};
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "finsymbden requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "finsymbden produces less output arguments");
/* ------------------------------------------------------------
   Get inputs LAD,perm,dz,firstq
   n = total number of dense columns/blocks, i.e. #cols in LAD
   m = number of constraints
   nperm = n - number of removed lorentz trace columns.
   ------------------------------------------------------------ */
  mxAssert(mxIsSparse(LAD_IN), "LAD must be sparse");           /* LAD */
  m = mxGetM(LAD_IN);
  n = mxGetN(LAD_IN);
  LADjc = mxGetJc(LAD_IN);
  LADir = mxGetIr(LAD_IN);
  permPr = mxGetPr(PERM_IN);                    /* perm */
  nperm = mxGetM(PERM_IN) * mxGetN(PERM_IN);
  dzJc = mxGetJc(DZ_IN);                        /* dz */
  dzIr = mxGetIr(DZ_IN);
  mxAssert(mxGetM(DZ_IN) == m && mxGetN(DZ_IN) == nperm, "dz size mismatch");
/* ------------------------------------------------------------
   INPUT firstQ == dense.l+1, points to 1st entry in dense.cols
    dealing with Lorentz-trace entries. Let lastQ point just beyond
    Lorentz trace/block entries, i.e. add n-nperm.
   ------------------------------------------------------------ */
  firstQ = (mwIndex) mxGetScalar(FIRSTQ_IN);         /*firstq, F-double to C-mwIndex.*/
  mxAssert(firstQ>0,"");
  --firstQ;
  lastQ = firstQ + n - nperm;
/* ------------------------------------------------------------
   Allocate integer working arrays:
   invdz(m), firstpiv(n), perm(n)
   ------------------------------------------------------------ */
  invdz = (mwIndex *) mxCalloc(MAX(1,m), sizeof(mwIndex));
  firstpiv = (mwIndex *) mxCalloc(MAX(1,n), sizeof(mwIndex));
  perm = (mwIndex *) mxCalloc(MAX(1,n), sizeof(mwIndex));
/* ------------------------------------------------------------
   Allocate OUTPUT mwIndex array dznewJc(n+1)
   ------------------------------------------------------------ */
  dznewJc = (mwIndex *) mxCalloc(n+1, sizeof(mwIndex));
/* ------------------------------------------------------------
   Let invdz(dzIr) = 1:nnz(dz). Note that nnz(dz)<m is the number
   subscripts that are actually in use.
   ------------------------------------------------------------ */
  nnzdz = dzJc[nperm];
  for(i = dzJc[0]; i < nnzdz; i++)
    invdz[dzIr[i]] = i;                 /* dz is m x nperm */
/* ------------------------------------------------------------
   Create new perm and dz-column pointers, to include lorentz trace cols.
   These cols are attached to Lorentz-blocks cols, whose subscripts
   range in firstQ:lastQ-1.
   ------------------------------------------------------------ */
  inz = 0;
  for(i = 0; i < nperm; i++){
    j = (mwIndex) permPr[i];
    perm[inz] = --j;
    dznewJc[inz++] = dzJc[i];
    if(j >= firstQ && j < lastQ){
/* ------------------------------------------------------------
   Attach Lorentz trace col. These cols are at nperm:n-1.
   ------------------------------------------------------------ */
      perm[inz] = nperm + j - firstQ;   /* insert associated trace column */
      mxAssert(perm[inz] < n,"");
      dznewJc[inz++] = dzJc[i+1];      /* no extra subscripts->start at end */
    }
  }
  mxAssert(inz == n,"");
  dznewJc[n] = dzJc[nperm];
/* ------------------------------------------------------------
   Compute firstpiv
   ------------------------------------------------------------ */
  getfirstpiv(firstpiv, invdz, dznewJc, LADjc,LADir, n);
/* ------------------------------------------------------------
   Outputs Lden.(LAD, perm, dz, first)
   ------------------------------------------------------------ */
  LDEN_OUT = mxCreateStructMatrix((mwSize)1,(mwSize)1, NLDEN_FIELDS, LdenFieldnames);
  LDEN_FIELD = mxDuplicateArray(LAD_IN);               /* LAD */
  mxSetField(LDEN_OUT,(mwIndex)0,"LAD",LDEN_FIELD);
  LDEN_FIELD = mxCreateDoubleMatrix(n, (mwSize)1, mxREAL);     /* perm */
  mxSetField(LDEN_OUT,(mwIndex)0,"perm",LDEN_FIELD);
  permPr = mxGetPr(LDEN_FIELD);
  for(i = 0; i < n; i++)
    permPr[i] = perm[i] + 1.0;                       /* C-mwIndex to F-double */
  LDEN_FIELD = mxDuplicateArray(DZ_IN);                /* dz */
/* NOTE: here we replace jc by dznewJc */
  mxFree(mxGetJc(LDEN_FIELD));
  mxSetJc(LDEN_FIELD, dznewJc);
  mxSetN(LDEN_FIELD, n);
  mxSetField(LDEN_OUT,(mwIndex)0,"dz",LDEN_FIELD);
  LDEN_FIELD  = mxCreateDoubleMatrix(n, (mwSize)1, mxREAL);  /* first */
  mxSetField(LDEN_OUT,(mwIndex)0,"first",LDEN_FIELD);
  firstPr = mxGetPr(LDEN_FIELD);
  for(i = 0; i < n; i++)
    firstPr[i] = firstpiv[i] + 1.0;               /* C-mwIndex to F-double */
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(perm);
  mxFree(firstpiv);
  mxFree(invdz);
}
