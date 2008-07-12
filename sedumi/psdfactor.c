/*
%                                          [ux,ispos] = psdfactor(x,K)
% PSDFACTOR  UX'*UX Cholesky factorization
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function [ux,ispos] = psdfactor(x,K)

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
#include <math.h>
#include "mex.h"
#include "blksdp.h"

#define UX_OUT myplhs[0]
#define ISPOS_OUT myplhs[1]
#define NPAROUT 2

#define X_IN prhs[0]
#define K_IN prhs[1]
#define NPARIN 2

/* ============================================================
   PSD: CHOLESKY FACTORIZATION
   ============================================================ */

/* ************************************************************
   PROCEDURE cholnopiv - U'*U factorization for nxn matrix,
     without pivoting.
   INPUT
     n - order of matrix to be factored
   UPDATED
     u - Full nxn. Input: Matrix to be factored.
       Output: Cholesky factor, X=triu(U)'*triu(U). 
       NOTE: tril(U,-1) is not affected by this function.
   RETURNS:
     0 = "success", 1 = "X is NOT positive definite"
   ************************************************************ */
mwIndex cholnopiv(double *u,const mwIndex n)
{
  mwIndex i,j;
  double uij,ujj;
  double *ui,*uj;

  if(n < 1)
    return 0;
/* ------------------------------------------------------------
   Solve the columns of U, for j=0:n-1
   ------------------------------------------------------------ */
  for(uj = u, j = 0; j < n; j++, uj+=n){
/* ------------------------------------------------------------
   Solve "uij" from the identity
   uii * uij = xij - u(1:i-1,i)'*u(1:i-1,j)
   ------------------------------------------------------------ */
    for(ui = u, i = 0, ujj = 0.0; i < j; i++, ui+=n){
      uij = (uj[i] = (uj[i] - realdot(ui,uj,i)) / ui[i]);
      ujj += SQR(uij);
    }
    ujj = uj[j] - ujj;
/* ------------------------------------------------------------
   By now, "ujj" should contain the final u(j,j)^2. Check whether
   it is positive. If not, then X was not p.d., thus fail.
   ------------------------------------------------------------ */
    if(ujj <= 0.0)
      return 1;               /* X is not positive definite */
    uj[j] = sqrt(ujj);
  }
  return 0;   /* success */
}

/* ************************************************************
   PROCEDURE prpicholnopiv  -  Computes triu matrix U s.t. U'*U = X.
     U, X in SeDuMi's complex format: X = [RE X, IM X]. No pivoting.
   INPUT:
     x - full 2*(n x n), should be Hermitian.
     n - order of u, x.
   UPDATED:
     u,upi - Both full nxn. Input: Matrix to be factored (real, imaginary).
       Output: Cholesky factor, X=U'*U with U=triu(u) + i*triu(upi,1). 
       NOTE: tril(u,-1) and tril(upi) are not affected by this function.
   RETURNS:
     0 = "success", 1 = "X is NOT positive definite"
   ************************************************************ */
mwIndex prpicholnopiv(double *u, double *upi, const mwIndex n)
{
  mwIndex i,j;
  double uii,uij,ujj;
  double *ui,*uipi, *uj, *ujpi;

  if(n < 1)
    return 0;
/* ------------------------------------------------------------
   Solve the columns of U, for j=0:n-1
   ------------------------------------------------------------ */
  for(uj=u, ujpi=upi, j = 0; j < n; j++, uj+=n, ujpi+=n){
    ujj = 0.0;
    for(ui=u, uipi=upi, i = 0; i < j; i++, ui+=n, uipi+=n){
/* ------------------------------------------------------------
   Solve "uij" from the identity
   uii * uij = xij - u(1:i-1,i)'*u(1:i-1,j)
   ------------------------------------------------------------ */
      uii = ui[i];
      uij = uj[i];                                     /* real part */
      uij -= realdot(ui,uj,i) + realdot(uipi,ujpi,i);
      uj[i] = (uij /= uii);
      ujj += SQR(uij);
      uij = ujpi[i];                                  /* imaginary part */
      uij += realdot(uipi,uj,i) - realdot(ui,ujpi,i);
      ujpi[i] = (uij /= uii);
      ujj += SQR(uij);
    }
    ujj = uj[j] - ujj;
/* ------------------------------------------------------------
   By now, "ujj" should contain the final u(j,j)^2. Check whether
   it is positive. If not, then X was not p.d., thus fail.
   ------------------------------------------------------------ */
    if(ujj <= 0.0)
      return 1;               /* X is not positive definite */
    uj[j] = sqrt(ujj);
  }
  return 0;
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   [ux,ispos] = psdfactor(x,K);
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  mxArray *myplhs[NPAROUT];
  coneK cK;
  mwIndex k,nk,nksqr, sdplen,sdpdim,lenfull, ispos;
  const double *x;
  double *ux;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "psdfactor requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "psdfactor produces less output arguments");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Compute statistics: sdpdim = rdim+hdim, sdplen = sum(K.s).
   ------------------------------------------------------------ */
  lenfull = cK.lpN +  cK.qDim + cK.rDim + cK.hDim;
  sdpdim = cK.rDim + cK.hDim;
  sdplen = cK.rLen + cK.hLen;
/* ------------------------------------------------------------
   Get input vector x, skip LP + Lorentz part
   ------------------------------------------------------------ */
  x = mxGetPr(X_IN);
  if(mxGetM(X_IN) * mxGetN(X_IN) == lenfull)
    x += cK.lpN + cK.qDim;
  else mxAssert(mxGetM(X_IN) * mxGetN(X_IN) == sdpdim, "x size mismatch.");
/* ------------------------------------------------------------
   Allocate output UX(sdpdim), ispos(1).
   ------------------------------------------------------------ */
  UX_OUT = mxCreateDoubleMatrix(sdpdim, (mwSize)1, mxREAL);
  ux = mxGetPr(UX_OUT);
  ISPOS_OUT = mxCreateDoubleMatrix((mwSize)1,(mwSize)1,mxREAL);
/* ------------------------------------------------------------
   PSD: Cholesky factorization.
   Initialize  ispos = 1 and ux = x.
   ------------------------------------------------------------ */
  ispos = 1;
  memcpy(ux, x, sdpdim * sizeof(double));       /* copy real + complex */
  for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
    nk = cK.sdpNL[k];
/* ------------------------------------------------------------
   Attempt Cholesky on block k. Returns 1 if fail (i.e. not psd).
   ------------------------------------------------------------ */
    if(cholnopiv(ux,nk)){
      ispos = 0;
      break;
    }
    triu2sym(ux,nk);
    ux += SQR(nk);
  }
/* ------------------------------------------------------------
   Complex Hermitian PSD Cholesky factorization, no pivoting.
   ------------------------------------------------------------ */
  if(ispos)
    for(; k < cK.sdpN; k++){                    /* complex Hermitian */
      nk = cK.sdpNL[k];
      nksqr = SQR(nk);
      if(prpicholnopiv(ux,ux+nksqr,nk)){
        ispos = 0;
        break;
      }
      triu2herm(ux,ux+nksqr,nk);
      ux += 2 * nksqr;
    }
/* ------------------------------------------------------------
   Return parameter ispos
   ------------------------------------------------------------ */
  *mxGetPr(ISPOS_OUT) = ispos;
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  k = MAX(nlhs, 1);
  memcpy(plhs,myplhs, k * sizeof(mxArray *));
  for(; k < NPAROUT; k++)
    mxDestroyArray(myplhs[k]);
}
