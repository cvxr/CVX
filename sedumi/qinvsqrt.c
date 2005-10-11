/*
   y = qinvsqrt(x,qdetx,K)
   Computes the Lorentz component of x^(-1/2).

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
#include "triuaux.h"
#include "blksdp.h"

#define Y_OUT plhs[0]
#define NPAROUT 1

#define X_IN prhs[0]
#define QDETX_IN prhs[1]
#define K_IN prhs[2]

#define NPARIN 3

/* ************************************************************
   PROCEDURE powminhalf - Computes y = x^{-1/2} for Lorentz,
     given rdetx = sqrt(det(x)).
   INPUT
     x - length n input vector.
     rdetx - scalar, qdetx = sqrt(det(x)).
     n - order
   OUTPUT
     y - length n vector, y = x^{-1/2}.
   ************************************************************ */
void powminhalf(double *y, const double *x, const double rdetx, const int n)
{
  double denom, alpha;
/* ------------------------------------------------------------
   Let denom = sqrt(sqrt(2)*x(1) + 2*rdetx), alpha = rdetx * denom
   ------------------------------------------------------------ */
  denom = sqrt(M_SQRT2 * x[0] + 2*rdetx);
  alpha = rdetx * denom;
/* ------------------------------------------------------------
   y(1) = x(1)/alpha + sqrt(2)/denom,
   y(2:n) = -x(2:n)/alpha
   ------------------------------------------------------------ */
  y[0] = x[0] / alpha + M_SQRT2 / denom;
  scalarmul(y+1, -1.0 / alpha, x+1, n-1);
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
  int i,k, nk;
  double *y;
  const double *x,*qdetx;
  coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "qinvsqrt requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "qinvsqrt generates less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get inputs x, qdetx
   ------------------------------------------------------------ */
  mxAssert(mxGetM(X_IN) * mxGetN(X_IN) >= cK.lpN + cK.qDim, "x size mismatch");
  x = mxGetPr(X_IN) + cK.lpN;              /* skip LP part */
  mxAssert(mxGetM(QDETX_IN) * mxGetN(QDETX_IN) == cK.lorN, "qdetx size mismatch");
  qdetx = mxGetPr(QDETX_IN);
/* ------------------------------------------------------------
   Allocate output y
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(cK.qDim, 1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   The actual job is done here: y = w^{-1/2}, Lorentz part.
   ------------------------------------------------------------ */
  for(k = 0; k < cK.lorN; k++){               /* LORENTZ */
    nk = cK.lorNL[k];
    powminhalf(y, x,qdetx[k],nk);
    y += nk; x += nk;
  }
}
