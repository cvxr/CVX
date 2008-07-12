/*
%                                                  y = invcholfac(u,K, perm)
% INVCHOLFAC  Computes y(perm,perm) = u' * u, with u upper triangular.
%
% SEE ALSO sedumi, getada3
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
function y = invcholfac(u,K, perm)

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
#include "mex.h"
#include "triuaux.h"

#define Y_OUT plhs[0]

#define NPAROUT 1

#define U_IN prhs[0]
#define K_IN prhs[1]
#define NPARINMIN 2
#define PERM_IN prhs[2]

#define NPARIN 3

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  mwIndex i, k, nk, nksqr, lenud, fwsiz;
  double *fwork, *fworkpi, *y, *permPr;
  const double *u;
  mwIndex *perm, *iwork;
  coneK cK;
  char isperm;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN){
    isperm = 0;
    mxAssert(nrhs >= NPARINMIN, "invcholfac requires at least 2 input arguments.");
  }
  else{
    permPr = mxGetPr(PERM_IN);
    isperm = (mxGetM(PERM_IN) * mxGetN(PERM_IN) > 0);
  }
  mxAssert(nlhs <= NPAROUT, "invcholfac generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  lenud = cK.rDim + cK.hDim;
/* ------------------------------------------------------------
   Get input U
   ------------------------------------------------------------ */
  u = mxGetPr(U_IN);
  mxAssert(mxGetM(U_IN) * mxGetN(U_IN) == lenud, "u size mismatch");
/* ------------------------------------------------------------
   Allocate output Y
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(lenud, (mwSize)1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   Allocate fwork = double( max(cK.rMaxn^2, 2*cK.hMaxn^2) )
   iwork = mwIndex(rLen+hLen)
   ------------------------------------------------------------ */
  fwsiz = MAX(SQR(cK.rMaxn),2*SQR(cK.hMaxn));
  fwork = (double *) mxCalloc( MAX(1,fwsiz), sizeof(double));
  iwork = (mwIndex *) mxCalloc( MAX(1, cK.rLen + cK.hLen), sizeof(mwIndex) );
/* ------------------------------------------------------------
   Let perm=iwork, and fworkpi = fwork + SQR(cK.hMaxn)
   ------------------------------------------------------------ */
  perm = iwork;
  fworkpi = fwork + SQR(cK.hMaxn);
/* ------------------------------------------------------------
   Convert Fortran to C-style in perm:
   ------------------------------------------------------------ */
  if(isperm){
    for(k = 0; k < cK.rLen + cK.hLen; k++){
      i = permPr[k];
      perm[k] = --i;
    }
/* ------------------------------------------------------------
   The actual job is done here:   Y = invperm(U'*U)
   ------------------------------------------------------------ */
    for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
      nk = cK.sdpNL[k];
      utmulx(fwork, u,u,nk);
      triu2sym(fwork,nk);
      invmatperm(y,fwork,perm,nk);                    /* Y(perm,perm) = Z */
      nksqr = SQR(nk);
      y += nksqr; u += nksqr;
      perm += nk;
    }
    for(; k < cK.sdpN; k++){
      nk = cK.sdpNL[k];
      nksqr = SQR(nk);
      prpiutmulx(fwork,fworkpi, u,u+nksqr,u,u+nksqr, nk);
      triu2herm(fwork,fworkpi,nk);
      invmatperm(y,fwork,perm,nk);                    /* Y(perm,perm) = Z */
      invmatperm(y+nksqr,fworkpi,perm,nk);        /* imaginary part */
      nksqr += nksqr;
      y += nksqr; u += nksqr;
      perm += nk;
    }
  }
  else{
/* ------------------------------------------------------------
   Without permutation, it's simply Y = U'*U.
   ------------------------------------------------------------ */
    for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
      nk = cK.sdpNL[k];
      utmulx(y, u,u,nk);
      triu2sym(y,nk);
      nksqr = SQR(nk);
      y += nksqr; u += nksqr;
    }
    for(; k < cK.sdpN; k++){
      nk = cK.sdpNL[k];
      nksqr = SQR(nk);
      prpiutmulx(y,y+nksqr, u,u+nksqr,u,u+nksqr, nk);
      triu2herm(y,y+nksqr,nk);
      nksqr += nksqr;
      y += nksqr; u += nksqr;
    }
  }
/* ------------------------------------------------------------
   Release working arrays.
   ------------------------------------------------------------ */
  mxFree(fwork);
  mxFree(iwork);
}
