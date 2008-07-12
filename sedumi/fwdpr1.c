/*
%                                                         y = fwdpr1(Lden, b)
% FWDPR1  Solves "PROD_k L(pk,betak) * y = b", where
%    where L(p,beta) = eye(n) + tril(p*beta',-1).
%
% SEE ALSO sedumi, dpr1fact, bwdpr1
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = fwdpr1(Lden, b)

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

/* y = fwdpr1(Lden, b) */
#define Y_OUT plhs[0]
#define NPAROUT 1

#define LDEN_IN prhs[0]
#define B_IN prhs[1]
#define NPARIN 2

/* ************************************************************
   PROCEDURE fwprodform - Solves PROD_j L(pj,betaj) * yNEW = yOLD.
   INPUT
     p - nonzeros of sparse m x n matrix P. Has xsuper(j+1) nonzeros in
      column j.
     xsuper - xsuper(j+1) is number of nonzeros in p(:,j).
     perm - lists pivot order for columns where ordered(j)==1.
     beta, betajc - Column beta(betajc[j]:betajc[j+1]-1) provides the
       beta vector in the j-th L-factor, L(p,beta)=eye(m) + tril(p*beta',-1).
     ordered - ordered[j]==1 iff p(:,j) and beta(L:,j) have been reordered;
       the original row numbers are in perm(:,j).
     n - number of columns in p, beta.
   UPDATED
     y - length m vector. On input, the rhs. On output the solution to
       PROD_j L(pj,betaj) * yNEW = yOLD.
   ************************************************************ */
void fwprodform(double *y, const mwIndex *xsuper, const mwIndex *perm,
                const double *p, const double *beta, const mwIndex *betajc,
                const char *ordered, const mwIndex n)
{
  mwIndex k,nk, mk;
/* ------------------------------------------------------------
   Forward solve L(pk,betak) * yNEXT = yPREV   for k=0,1,..,n-1, resp.
   ------------------------------------------------------------ */
  for(k = 1; k <= n; k++){
    mk = xsuper[k];
    nk = betajc[k] - betajc[k-1];
    if(ordered[k-1]){
      fwipr1o(y, perm, p, beta, mk, nk);
      perm += mk;
    }
    else
      fwipr1(y, p, beta, mk, nk);
    beta += nk;
    p += mk;
  }
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
       y = fwdpr1(Lden, b)
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
 const mxArray *MY_FIELD;
 char *ordered;
 mwIndex m,n,nden,dznnz, i,j, permnnz;
 const double *beta, *betajcPr, *orderedPr, *pivpermPr, *p;
 mwIndex *betajc, *pivperm;
 double *y, *fwork;
 jcir dz;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "fwdpr1 requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "fwdpr1 generates less output arguments.");
/* ------------------------------------------------------------
   Get input b
   ------------------------------------------------------------ */
  mxAssert(!mxIsSparse(B_IN), "b should be full");
  m = mxGetM(B_IN);
  n = mxGetN(B_IN);
/* ------------------------------------------------------------
   Create output y as a copy of b
   ------------------------------------------------------------ */
  Y_OUT = mxDuplicateArray(B_IN);                 /* Y_OUT = B_IN */
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   Disassemble dense-update structure Lden (p,xsuper,beta,betajc,rowperm)
   NOTE: if there are no dense columns, then simply let y = b, and return.
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(LDEN_IN), "Parameter `Lden' should be a structure.");
  MY_FIELD = mxGetField(LDEN_IN,(mwIndex)0,"betajc");
  mxAssert( MY_FIELD!= NULL, "Missing field Lden.betajc.");    /* betajc */
  nden = mxGetM(MY_FIELD) * mxGetN(MY_FIELD) - 1;
/* If no dense columns: get out immediately ! */
  if(nden == 0)
    return;
  else{
    betajcPr = mxGetPr(MY_FIELD);
    MY_FIELD = mxGetField(LDEN_IN,(mwIndex)0,"p");
    mxAssert( MY_FIELD != NULL, "Missing field Lden.p.");         /* p */
    p = mxGetPr(MY_FIELD);
    MY_FIELD = mxGetField(LDEN_IN,(mwIndex)0,"dopiv");
    mxAssert( MY_FIELD != NULL, "Missing field Lden.dopiv.");  /* dopiv */
    mxAssert(mxGetM(MY_FIELD) * mxGetN(MY_FIELD) == nden, "Size mismatch Lden.dopiv.");
    orderedPr = mxGetPr(MY_FIELD);
    MY_FIELD = mxGetField(LDEN_IN,(mwIndex)0,"pivperm"); /* pivperm */
    mxAssert( MY_FIELD != NULL, "Missing field Lden.pivperm.");
    pivpermPr = mxGetPr(MY_FIELD);
    permnnz = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
    MY_FIELD = mxGetField(LDEN_IN,(mwIndex)0,"beta");    /* beta */
    mxAssert( MY_FIELD != NULL, "Missing field Lden.beta.");
    beta = mxGetPr(MY_FIELD);
    MY_FIELD = mxGetField(LDEN_IN,(mwIndex)0,"dz");      /* dz */
    mxAssert( MY_FIELD != NULL, "Missing field Lden.dz.");
    mxAssert(mxGetM(MY_FIELD) == m && mxGetN(MY_FIELD) == nden, "Lden.dz size mismatch.");
    mxAssert(mxIsSparse(MY_FIELD), "Lden.dz must be sparse.");
    dz.jc = mxGetJc(MY_FIELD);
    dz.ir = mxGetIr(MY_FIELD);                             /* (rowperm) */
    dznnz = dz.jc[nden];
    mxAssert(dznnz <= m, "Lden.dz size mismatch.");
/* ------------------------------------------------------------
   Allocate working arrays mwIndex: betajc(nden+1), pivperm(permnnz),
   char: ordered(nden)
   double: fwork(dznnz).
   ------------------------------------------------------------ */
    betajc = (mwIndex *) mxCalloc(nden + 1,sizeof(mwIndex));
    pivperm = (mwIndex *) mxCalloc(MAX(permnnz,1),sizeof(mwIndex));
    ordered = (char *) mxCalloc(nden,sizeof(char));   /* nden > 0 */
    fwork = (double *) mxCalloc(MAX(dznnz,1),sizeof(double));
/* ------------------------------------------------------------
   Convert betajcPr, ordered, pivperm to mwIndex
   ------------------------------------------------------------ */
    for(i = 0; i <= nden; i++){
      j = betajcPr[i];
      betajc[i] = --j;
    }
    for(i = 0; i < nden; i++){
      ordered[i] = orderedPr[i];
    }
    for(i = 0; i < permnnz; i++){
      pivperm[i] = pivpermPr[i];
    }
    MY_FIELD = mxGetField(LDEN_IN,(mwIndex)0,"beta");
    mxAssert(mxGetM(MY_FIELD) * mxGetN(MY_FIELD) == betajc[nden], "Size mismatch Lden.beta.");
/* ------------------------------------------------------------
   The actual job is done here: fwork = PROD_L\b(perm,j), y(perm,j) = fwork.
   ------------------------------------------------------------ */
    for(j = 0; j < n; j++, y += m){
      for(i = 0; i < dznnz; i++)            /* fwork = y(dzir) */
        fwork[i] = y[dz.ir[i]];
      fwprodform(fwork, dz.jc, pivperm, p, beta, betajc, ordered, nden);
      for(i = 0; i < dznnz; i++)            /* y(dzir) = fwork */
        y[dz.ir[i]] = fwork[i];
    }
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
    mxFree(fwork);
    mxFree(ordered);
    mxFree(pivperm);
    mxFree(betajc);
  } /* if dense columns */
}
