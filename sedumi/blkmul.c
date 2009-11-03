/*
    y = blkmul(mu,d,nL)
    yields length(y)=sum(nL) vector with y[k] = mu(k) * d[k].
    length(mu)=length(nL), length(d)=sum(nL).

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
#include "blksdp.h"

#define Y_OUT plhs[0]
#define NPAROUT 1
#define MU_IN prhs[0]
#define D_IN prhs[1]
#define NL_IN prhs[2]
#define NPARIN 3

/* ************************************************************
   PROCEDURE blkmul - Multiply each block "dk" of length nL(k)
     with scalar mu(k), yielding "yk".

   INPUT
     mu kappa-vector
     d n-vector
     nL kappa-array, with sum(nL)=n. length of each block k=1:kappa.
     kappa, n - order of mu, d resp.
   OUTPUT
     y n-vector. On output, y(1:nL(1))= kappa(1)*d(1:nL(1)), etc.
   ************************************************************ */
int blkmul(double *y, const double *mu,const double *d,const mwIndex *nL,
           const mwIndex kappa, mwIndex n)
{
  mwIndex k,nk;
/* ------------------------------------------------------------
   For each block k: yk = mu(k) * d[k].
   n denotes remaining length of d and y.
   ------------------------------------------------------------ */
  for(k = 0; k < kappa; k++){
    nk = nL[k];
    if(n < nk)
      return 1;
    scalarmul(y,mu[k],d,nk);
    y += nk;
    d += nk;
    n -= nk;
  }
/* ------------------------------------------------------------
   If n = sum(nL) then now n=0
   ------------------------------------------------------------ */
  return n;
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
  mwIndex i,kappa,n;
  double *y;
  const double *d, *mu, *nLpr;
  mwIndex *nL;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
 mxAssert(nrhs >= NPARIN, "blkmul requires more input arguments.");
 mxAssert(nlhs <= NPAROUT, "blkmul generates 1 output argument.");
/* ------------------------------------------------------------
   Get inputs d, mu, NL
   ------------------------------------------------------------ */
 d = mxGetPr(D_IN);
 n = mxGetM(D_IN) * mxGetN(D_IN);
 mu = mxGetPr(MU_IN);
 kappa = mxGetM(MU_IN) * mxGetN(MU_IN);
 mxAssert(mxGetM(NL_IN) * mxGetN(NL_IN) == kappa, "nL size mismatch");
 nLpr = mxGetPr(NL_IN);
/* ------------------------------------------------------------
   Allocate output y(qDim)
   ------------------------------------------------------------ */
 Y_OUT =  mxCreateDoubleMatrix(n, 1, mxREAL);
 y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   ALLOCATE working array nL(kappa) = nLPr in ints.
   ------------------------------------------------------------ */
 nL = (mwindex *) mxCalloc(MAX(1,kappa), sizeof(mwIndex));
 for(i = 0; i < kappa; i++)
   nL[i] = nLpr[i];
/* ------------------------------------------------------------
   Multiply each block "dk" of length nL(k) with scalar mu(k).
   ------------------------------------------------------------ */
 blkmulsize=blkmul(y, mu,d,nL,kappa,n);
 mxAssert(blkmulsize == 0, "nL size mismatch");
/* ------------------------------------------------------------
   RELEASE working array nL
   ------------------------------------------------------------ */
 mxFree(nL);
}
