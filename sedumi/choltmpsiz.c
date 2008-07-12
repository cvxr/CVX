/*
 L.tmpsiz = choltmpsiz(L)

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

#define TMPSIZ_OUT plhs[0]
#define NPAROUT 1

#define L_IN prhs[0]       /* symbolic Cholesky structure: {L.L, L.xsuper} */
#define NPARIN 1

#include "mex.h"

/* ************************************************************
   PROCEDURE gettmpsiz - Compute "fwork"-size in PRECORRECT of
     BLKCHOL2. Since fwork = -(Lk * inv(Dk) * Lk')_{snode j}, we have
     tmpsiz = MAX_{k,j in SUPER} mk * q - q(q-1)/2,
     with q := ncolup(k,j) = #nz-rows in L(:,k) corresponding to
               subnodes of j.
          mk := #nz-rows in L(:,k) corresponding to subnodes of j : nsuper.
   INPUT
     ljc,lir - sparsity structure of m x m matrix L (not compressed).
     xsuper,nsuper - supernodal partition of nodes 1:m.
     snode - maps nodes 1:m to supernode containing it.
   RETURNS tmpsiz.
   ************************************************************ */
mwIndex gettmpsiz(const mwIndex *ljc,const mwIndex *lir,const mwIndex *xsuper,
              const mwIndex nsuper, const mwIndex *snode)
{
  mwIndex tmpsiz, ksup,k,j,i,nextj, mk,inz,ncolup, sizkj,ubsiz;

  tmpsiz = 0;
/* ------------------------------------------------------------
   For each supernode ksup = 1:nsuper, and affected supernode snode[j]:
   ncolup = #nz-rows in L(:,k) corresponding to subnodes of snode[j].
   mk := #nz-rows in L(:,k) corresponding to subnodes of snode[j] : nsuper.
   ------------------------------------------------------------ */
  for(ksup = 0; ksup < nsuper; ksup++){
    k = xsuper[ksup];
/* ------------------------------------------------------------
   Let mk be number of below-diag-block(k) nonzeros. This
   is upper bound on both q and mk.
   ------------------------------------------------------------ */
    inz = ljc[k] + (xsuper[ksup+1] - k);        /* points below diag-block */
    mk = ljc[k+1] - inz;
    ubsiz = mk * (mk + 1) / 2;                  /* ubound on tmpsiz(k) */
    i = lir[ljc[k+1]-1];                        /* last subscript in k */
/* ------------------------------------------------------------
   Browse through all affected supernodes snode[j], as long as
   they're worth considering (i.e. ubsiz > tmpsiz).
   ------------------------------------------------------------ */
    while((inz < ljc[k+1]) && (ubsiz > tmpsiz)){
      j = lir[inz];                           /* 1st affected column */
      nextj = xsuper[snode[j] + 1];           /* beyond supernode of j */
      if(i < nextj){
        ncolup = mk;                          /* last affected supernode */
      }
      else{
        ncolup = 1;                        /* Compute #affected cols in j */
        for(++inz; lir[inz] < nextj; inz++)
          ncolup++;
      }
      sizkj = mk * ncolup - ncolup*(ncolup-1)/2;
      if(sizkj > tmpsiz)
        tmpsiz = sizkj;
      mk -= ncolup;                       /* proceed beyond snode[j] */
      ubsiz = mk * (mk + 1) / 2;          /* ubound on tmpsiz(k,[j+1,:]) */
    } /* inz in column L(:,k) AND ubsiz > tmpsiz */
  } /* for all supernodes ksup */
  return tmpsiz;
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
  mwIndex i,j, nsuper,m, jsup,tmpsiz;
  const mwIndex *ljc,*lir;
  mwIndex *xsuper, *snode;
  const double *xsuperPr;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "choltmpsiz requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "choltmpsiz produces less output arguments.");
/* ------------------------------------------------------------
   Disassemble block Cholesky structure L
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(L_IN), "Parameter `L' should be a structure.");
  L_FIELD = mxGetField(L_IN,(mwIndex)0,"L");          /* L.L */
  mxAssert( L_FIELD != NULL, "Missing field L.L.");
  m = mxGetM(L_FIELD);
  mxAssert(m == mxGetN(L_FIELD), "L.L must be square.");
  mxAssert(mxIsSparse(L_FIELD), "L.L should be sparse.");
  ljc = mxGetJc(L_FIELD);
  lir = mxGetIr(L_FIELD);
  L_FIELD = mxGetField(L_IN,(mwIndex)0,"xsuper");         /* L.xsuper */
  mxAssert( L_FIELD != NULL, "Missing field L.xsuper.");
  nsuper = mxGetM(L_FIELD) * mxGetN(L_FIELD) - 1;
  mxAssert( nsuper <= m, "Size L.xsuper mismatch.");
  xsuperPr = mxGetPr(L_FIELD);
/* ------------------------------------------------------------
   Allocate working arrays:
   ------------------------------------------------------------ */
  xsuper    = (mwIndex *) mxCalloc(nsuper+1,sizeof(mwIndex));
  snode     = (mwIndex *) mxCalloc(m,sizeof(mwIndex));
/* ------------------------------------------------------------
   Convert XSUPER to integer and C-Style
   ------------------------------------------------------------ */
  for(i = 0; i <= nsuper; i++){
    j =  (mwIndex) xsuperPr[i];
    mxAssert(j>0,"");
    xsuper[i] = --j;
  }
/* ------------------------------------------------------------
   SNODE: map each column to the supernode containing it
   ------------------------------------------------------------ */
  j = xsuper[0];
  for(jsup = 0; jsup < nsuper; jsup++){
    while(j < xsuper[jsup + 1])
      snode[j++] = jsup;
  }
/* ------------------------------------------------------------
   The main job: compute (upper bound on) blkchol-tmpsiz.
   ------------------------------------------------------------ */
  tmpsiz = gettmpsiz(ljc,lir,xsuper,nsuper, snode);
/* ------------------------------------------------------------
   return OUTPUT variable tmpsiz
   ------------------------------------------------------------ */
  TMPSIZ_OUT = mxCreateDoubleMatrix((mwSize)1,(mwSize)1,mxREAL);          /* L.tmpsiz */
  *mxGetPr(TMPSIZ_OUT) = (double) tmpsiz;
/* ------------------------------------------------------------
   Release working arrays.
   ------------------------------------------------------------ */
  mxFree(snode);
  mxFree(xsuper);
}
