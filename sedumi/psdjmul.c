/*
%                                                    y = psdjmul(x,y, K)
% PSDJMUL  for full x,y. Computes (XY+YX)/2
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = psdjmul(x,y, K)

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

#define Z_OUT plhs[0]
#define NPAROUT 1

#define X_IN prhs[0]
#define Y_IN prhs[1]
#define K_IN prhs[2]
#define NPARIN 3

/* ============================================================ 
   PSD: Z = tril(X*Y+Y*X)/2
   ============================================================ */

/* ************************************************************
   PROCEDURE symjmul(z,x,y,n)  --
     Z = tril(X * Y + Y * X) / 2, with X,Y symmetric.
     The strict upper triangular of Z is undefined.
   INPUT
     x,y - symmetric n x n
     n - order of square x,y,z matrices
   UPDATED
     z - full n x n, on entry arbitrary, on return tril(Z) = tril(XY+YX)/2.
   ************************************************************ */
void symjmul(double *z,const double *x,const double *y,const mwIndex n)
{
 mwIndex i,j,icol;

 /* ------------------------------------------------------------
    For i=0..n-1:
        for j=i..n-1:  let z(j,i) = (x(:,i)'*y(:,j) + y(:,i)'*x(:,j))/2
    ------------------------------------------------------------ */
 for(j = 0; j < n; z += n, x+=n, y+=n, j++){
   z[j] = realdot(x,y,n);
   for(i = j+1, icol = n; i < n; icol+=n, i++)
     z[i] = (realdot(x+icol,y,n) + realdot(x,y+icol,n)) / 2;
 }
}

/* ************************************************************
   PROCEDURES hermjmul --
     Z = tril(X * Y + Y * X) / 2, with X,Y Hermitian. X = [RE X,IM X].
     Only stict lower triangular and real diagonal of R are defined.
   INPUT
     x,y - Hermitian; full 2*(n x n).
     n - order of square x,y,z matrices
   UPDATED
     z - full 2*(n x n), on entry arbitrary, on return
       tril(Z) = tril(XY + YX)/2.
   ************************************************************ */
void hermjmul(double *z,const double *x,const double *y,const mwIndex n)
{
 mwIndex i,j,icol,jcol,nsqr;
 const double *xpi,*ypi;

 nsqr = SQR(n);
 xpi = x + nsqr; ypi = y + nsqr;
 /* ------------------------------------------------------------
    For i=0..n-1:
        for j=i..n-1:  let z(j,i) = (x(:,i)'*y(:,j) + y(:,i)'*x(:,j))/2
    ------------------------------------------------------------ */
 for(j = 0, jcol=0; j < n; z += n, jcol+=n, j++){
   z[j] = realdot(x+jcol,y+jcol,n) + realdot(xpi+jcol,ypi+jcol,n);
   for(i = j+1, icol = jcol; i < n; i++){
     icol += n;
     z[i] = (realdot(x+icol,y+jcol,n) + realdot(x+jcol,y+icol,n)
       + realdot(xpi+icol,ypi+jcol,n) + realdot(xpi+jcol,ypi+icol,n)) / 2;
   }
 }
 for(j = 0, jcol=0; j < n; z += n, jcol+=n, j++){
   for(i = j+1, icol = jcol; i < n; i++){
     icol += n;
     z[i] = (realdot(x+icol,ypi+jcol,n)
	     - realdot(xpi+icol,y+jcol,n)
	     - realdot(x+jcol,ypi+icol,n)
	     + realdot(xpi+jcol,y+icol,n)) / 2;
   }
 }
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
    z = psdjmul(x,y, K)
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
 mwIndex k, nk, nksqr, lenfull, ifirst, lenud;
 double *z;
 const double *x,*y;
 coneK cK;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
 mxAssert(nrhs >= NPARIN, "jmulK requires more input arguments.");
 mxAssert(nlhs <= NPAROUT, "jmulK generates 1 output argument.");
 /* ------------------------------------------------------------
    Disassemble cone K structure
    ------------------------------------------------------------ */
 conepars(K_IN, &cK);
 /* ------------------------------------------------------------
    Get statistics of cone K structure
    ------------------------------------------------------------ */
 ifirst = cK.lpN + cK.qDim;          /* point to PSD */
 lenud = cK.rDim + cK.hDim;
 lenfull = ifirst + lenud;
/* ------------------------------------------------------------
   Get inputs x and y.
   ------------------------------------------------------------ */
 x = mxGetPr(X_IN) + ifirst;
 y = mxGetPr(Y_IN) + ifirst;
 mxAssert(!mxIsSparse(X_IN) && !mxIsSparse(Y_IN), "Sparse inputs not supported by this version of jmulK.");
 mxAssert(mxGetM(X_IN) * mxGetN(X_IN) == lenfull, "size x mismatch.");
 mxAssert(mxGetM(Y_IN) * mxGetN(Y_IN) == lenfull, "size y mismatch.");
/* ------------------------------------------------------------
   Allocate output Z
   ------------------------------------------------------------ */
 Z_OUT =  mxCreateDoubleMatrix(lenud, (mwIndex)1, mxREAL);
 z = mxGetPr(Z_OUT);
/* ------------------------------------------------------------
   The actual job is done here: Z = (XY + YX)/2
   ------------------------------------------------------------ */
 for(k=0; k < cK.rsdpN; k++){                /* real symmetric */
   nk = cK.sdpNL[k];                 /* real to mwIndex */
   symjmul(z,x,y,nk);
   tril2sym(z,nk);
   nksqr = SQR(nk);
   z += nksqr; x += nksqr;
   y += nksqr;
 }
 for(; k < cK.sdpN; k++){                    /* complex Hermitian */
   nk = cK.sdpNL[k];                 /* real to mwIndex */
   nksqr = SQR(nk);
   hermjmul(z,x,y,nk);
   tril2herm(z,z+nksqr,nk);
   nksqr += nksqr;
   z += nksqr; x += nksqr;
   y += nksqr;
 }
}
