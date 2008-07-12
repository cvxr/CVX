/*
%                                                 y = triumtriu(r,u,K)
% TRIUMTRIU  Computes y = r * u
%   Both r and u should be upper triangular.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = triumtriu(r,u,K)

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
#include "blksdp.h"

#define Y_OUT plhs[0]

#define NPAROUT 1

#define R_IN prhs[0]
#define U_IN prhs[1]
#define K_IN prhs[2]

#define NPARIN 3

/* ************************************************************
   PROCEDURE realumulu - Computes Y = R*U with R, U upper triu.
   Y will also be triu.
   ************************************************************ */
void realumulu(double *y, const double *r,const double *u,const mwIndex n)
{
  mwIndex j,jcol,k;

/* ------------------------------------------------------------
   For each column Rk, we have Y += Rk * U(k,:) = R(1:k,k)*U(k,k:n).
   Let y,r point to y(1,k), r(1,k), and u to u(k,k).
   ------------------------------------------------------------ */
  for(k = 0; k < n; k++, y += n, r += n, u += n+1)    /* each column Rk: */
    for(j = k, jcol = 0; j < n; j++, jcol += n)    /* each j >= k: */
      addscalarmul(y+jcol, u[jcol], r, k+1);       /* y(:,j) += r(:,k) * ukj */
}

/* ************************************************************
   PROCEDURE prpiumulu - Computes Y = R*U with R, U upper triu,
    real diagonal. Y will also be triu, real diag.
   ************************************************************ */
void prpiumulu(double *y,double *ypi, const double *r,const double *rpi,
               const double *u,const double *upi, const mwIndex n)
{
  mwIndex j,jcol,k;

/* ------------------------------------------------------------
   For each column Rk, we have Y += Rk * U(k,:) = R(1:k,k)*U(k,k:n).
   Let y,r point to y(1,k), r(1,k), and u to u(k,k).
   ------------------------------------------------------------ */
  for(k = 0; k < n; k++, y+=n,ypi+=n, r+=n,rpi+=n, u+=n+1,upi+=n+1){
/* y(:,k) += Rk * ukk, with ukk and rkk real. */
    addscalarmul(y, u[0], r, k+1);
    addscalarmul(ypi, u[0], rpi, k);
    for(j = k+1, jcol = n; j < n; j++, jcol += n){    /* each j > k: */
/* y(:,j) += Rk * ukj, with rkk real. */
      addscalarmul(y+jcol, u[jcol], r, k+1);          /* real part */
      addscalarmul(y+jcol, -upi[jcol], rpi, k);
      addscalarmul(ypi+jcol, u[jcol], rpi, k);        /* imaginary part */
      addscalarmul(ypi+jcol, upi[jcol], r, k+1);
    }
  }
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
  mwIndex k, nk, nksqr, lenud;
  double *y;
  const double *r,*u;
  coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "triumtriu requires at least 3 input arguments.");
  mxAssert(nlhs <= NPAROUT, "triumtriu generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  lenud = cK.rDim + cK.hDim;
/* ------------------------------------------------------------
   Get input R, U
   ------------------------------------------------------------ */
  mxAssert(mxGetM(R_IN) * mxGetN(R_IN) == lenud, "r size mismatch");
  mxAssert(mxGetM(U_IN) * mxGetN(U_IN) == lenud, "u size mismatch");
  r = mxGetPr(R_IN);
  u = mxGetPr(U_IN);
/* ------------------------------------------------------------
   Allocate output Y
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(lenud, (mwSize)1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   The actual job is done here: Y = R*U  (triu * triu)
   ------------------------------------------------------------ */
  for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
    nk = cK.sdpNL[k];
    realumulu(y, r,u,nk);
    triu2sym(y,nk);                            /* Give also tril-factor */
    nksqr = SQR(nk);
    y += nksqr; r += nksqr; u += nksqr;
  }
  for(; k < cK.sdpN; k++){                    /* complex Hermitian */
    nk = cK.sdpNL[k];
    nksqr = SQR(nk);
    prpiumulu(y,y+nksqr, r,r+nksqr, u,u+nksqr, nk);
    triu2herm(y,y+nksqr,nk);                 /* Give also tril-factor */
    nksqr += nksqr;
    y += nksqr; r += nksqr; u += nksqr;
  }
}
