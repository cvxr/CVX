/*
%                                          z = psdinvjmul(xlab,xfrm, y, K)
% PSDINVJMUL  solves x jmul z = y, with x = XFRM*diag(xlab)*XFRM'
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function z = psdinvjmul(xlab,xfrm, y, K)

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
#include "reflect.h"

/* z = psdinvjmul(xlab,xfrm, y, K) */
#define Z_OUT plhs[0]
#define NPAROUT 1

#define X_IN prhs[0]
#define FRM_IN prhs[1]
#define Y_IN prhs[2]
#define K_IN prhs[3]
#define NPARIN 4

/* ============================================================ 
   PSD: z(i.j) = 2*y(i,j)/(xi + xj)   (Z and Y are matrices, x is vector)
   ============================================================ */

/* ************************************************************
   PROCEDURE diagjdiv  -   yij := 2*y(i,j)/(xi + xj) for i >= j.
     The strict upper triangular of Y is left unchanged.
   INPUT
     x - length n vector
     n - order of x and square y matrix
   UPDATED
     y - full n x n, on return tril(yNEW) = 2*y(i,j)/(xi + xj), i >= j.
   ************************************************************ */
void diagjdiv(double *y,const double *x,const mwIndex n)
{
  mwIndex i,j;
  double xj;

/* ------------------------------------------------------------
   For j=0..n-1:
   for i=j..n-1:  let y(i,j) *= 2/(xi + xj)
   ------------------------------------------------------------ */
 for(j = 0; j < n; y+=n, j++){
   xj = x[j];
   y[j] /= xj;
   for(i = j+1; i < n; i++)
     y[i] *= 2 / (x[i] + xj);
 }
}

/* ************************************************************
   PROCEDURE psdinvjmul
   INPUT
     frms lenud+hLen-vector: contains coded orthogonal matrix "Q" for
        each PSD block, listed in sdpNL. Eigvecs of X.
     x sum(K.s)-vector: contains eigenvalues, so that X = Q*LAB*Q'.
     y lenud-vector: full matrix for each PSD block.
     sdpNL sdpN-array: order of each PSD block, i.e. K.s.
     rsdpN number of real symmetric blocks: 0<= rsdpN <= sdpN.
     sdpN number of PSD blocks, i.e. length(K.s).
   OUTPUT:
     z lenud-vector. On ouput, X*Z+Z*X = 2 * Y.
   WORKING ARRAY:
     fwork length max(rmaxn,2*hmaxn).
   ************************************************************ */
void psdinvjmul(double *z, const double *frms, const double *x,
                const double *y, const mwIndex *sdpNL,const mwIndex rsdpN,
                const mwIndex sdpN, double *fwork)
{
  mwIndex k,nk,nksqr;
  const double *beta;
/* ------------------------------------------------------------
   PSD: Since X = Q'*diag(x)*Q, we have XZ+ZX = 2Y iff
    qzqt(i,j) = 2*qyqt(i,j)/(xi + xj),  qzqt = Q*z*Q', qyqt = Q*y*Q'.
   ------------------------------------------------------------ */
  for(k = 0; k < rsdpN; k++){
   nk = sdpNL[k];
   nksqr = SQR(nk);
/* ------------------------------------------------------------
   Let z = Q*y*Q'
   ------------------------------------------------------------ */
   memcpy(z, y, nksqr * sizeof(double));
   beta = frms + nksqr - nk;               /* beta = frms(:,end) */
   qxqt(z, beta, frms, nk, fwork);

/* ------------------------------------------------------------
   Solve diag(x) jmul zNEW = zOLD
   ------------------------------------------------------------ */
   diagjdiv(z,x,nk);
/* ------------------------------------------------------------
   Let zFINAL = Q'*z*Q  (back into old eig-basis)
   ------------------------------------------------------------ */
   qtxq(z, beta, frms, nk, fwork);
   tril2sym(z,nk);
   z += nksqr; y += nksqr; frms += nksqr;
   x += nk;
 }
 for(; k < sdpN; k++){                    /* complex Hermitian */
   nk = sdpNL[k];
   nksqr = SQR(nk);
/* ------------------------------------------------------------
   Let z = Q*y*Q'
   ------------------------------------------------------------ */
   memcpy(z, y, 2 * nksqr * sizeof(double));
   beta = frms + 2 * nksqr;                   /* beta = frms(:,2*n+1) */
   prpiqxqt(z,z+nksqr, beta, frms,frms+nksqr, nk, fwork);
/* ------------------------------------------------------------
   Solve diag(x) jmul zNEW = zOLD. Since x is real, we can handle
   the real and imaginary parts seperately.
   ------------------------------------------------------------ */
   diagjdiv(z,x,nk);            /* real part */
   diagjdiv(z+nksqr,x,nk);      /* imaginary part */
/* ------------------------------------------------------------
   Let zFINAL = Q'*z*Q  (back into old eig-basis)
   ------------------------------------------------------------ */
   prpiqtxq(z,z+nksqr, beta, frms,frms+nksqr, nk, fwork);
   tril2herm(z,z+nksqr,nk);
   nksqr += nksqr;
   z += nksqr; y += nksqr; frms += nksqr + nk;          /* skip also beta. */
   x += nk;
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
 mwIndex i, lenfull, lendiag, lenud, qsize;
 double *z, *fwork;
 const double *x,*y, *frms;
 mwIndex *sdpNL;
 coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
 mxAssert(nrhs >= NPARIN, "psdinvjmul requires more input arguments.");
 mxAssert(nlhs <= NPAROUT, "psdinvjmul generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
 conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
 lenud = cK.rDim + cK.hDim;
 qsize = lenud + cK.hLen;
 lenfull = cK.lpN +  cK.qDim + lenud;
 lendiag = cK.lpN + 2 * cK.lorN + cK.rLen + cK.hLen;
/* ------------------------------------------------------------
   Get inputs x, frm, y.
   ------------------------------------------------------------ */
 mxAssert(!mxIsSparse(X_IN) && !mxIsSparse(Y_IN), "Sparse inputs not supported by this version of psdinvjmul.");
 x = mxGetPr(X_IN);                                /* get x and y */
 y = mxGetPr(Y_IN);
 if(mxGetM(Y_IN) * mxGetN(Y_IN) != lenud){
   mxAssert(mxGetM(Y_IN) * mxGetN(Y_IN) == lenfull, "size y mismatch.");
   y += cK.lpN + cK.qDim;          /* point to PSD */
 }
 if(mxGetM(X_IN) * mxGetN(X_IN) != cK.rLen + cK.hLen){
   mxAssert(mxGetM(X_IN) * mxGetN(X_IN) == lendiag, "size xlab mismatch.");
   x += cK.lpN + 2 * cK.lorN;       /* point to PSD */
 }
 mxAssert(mxGetM(FRM_IN) * mxGetN(FRM_IN) == qsize, "size xfrm mismatch.");
 frms = mxGetPr(FRM_IN);
/* ------------------------------------------------------------
   Allocate output Z
   ------------------------------------------------------------ */
 Z_OUT =  mxCreateDoubleMatrix(lenud, (mwIndex)1, mxREAL);
 z = mxGetPr(Z_OUT);
/* ------------------------------------------------------------
   Allocate working array fwork(max(rmaxn,2*hmaxn))
   integer working array sdpNL(sdpN).
   ------------------------------------------------------------ */
  fwork = (double *) mxCalloc(MAX(1,MAX(cK.rMaxn,2*cK.hMaxn)),sizeof(double));
  sdpNL = (mwIndex *) mxCalloc(MAX(1,cK.sdpN), sizeof(mwIndex));
/* ------------------------------------------------------------
   double to integer
   ------------------------------------------------------------ */
  for(i = 0; i < cK.sdpN; i++)
    sdpNL[i] = (mwIndex) cK.sdpNL[i];
  psdinvjmul(z,frms,x,y,sdpNL,cK.rsdpN,cK.sdpN, fwork);
/* ------------------------------------------------------------
   Release working array
   ------------------------------------------------------------ */
  mxFree(sdpNL);
  mxFree(fwork);
}
