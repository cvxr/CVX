/*
%                                                [lab,q] = psdeig(x,K)
% PSDEIG  Computes spectral coefficients of x w.r.t. K
%   Arguments "q" is optional - without it's considerably faster.
%   FLOPS indication: 1.3 nk^3 versus 9.0 nk^3 for nk=500,
%                     1.5 nk^3        9.8 nk^3 for nk=50.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function [lab,q] = psdeig(x,K)

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

#define USE_SVD 1 /* Uses [q,lab,q]=svd(X) for complex hermitian (psd?) x */

#define LAB_OUT plhs[0]
#define Q_OUT plhs[1]
#define NPAROUT 2

#define X_IN prhs[0]
#define K_IN prhs[1]
#define NPARIN 2

/* ============================================================
   PSD: projection onto symmetric/ skew-symmetric routines.
   ============================================================ */
/* ************************************************************
   PROCEDURE symproj -- Y = (X+X')/2
   INPUT x, n - full n x n matrix x.
   OUTPUT y - on output, contains (x+x')/2
   ************************************************************ */
void symproj(double *y, const double *x, const mwIndex n)
{
  mwIndex colp,i,j;
  double yij;

  /* ------------------------------------------------------------
     x points to x(:,i);     x+colp = x(:,j).
     ------------------------------------------------------------ */
  for(i = 0; i < n; x += n, y += n, i++){
    y[i] = x[i];                         /* diagonal entry */
    for(colp = n + i, j=i+1; j<n; j++, colp += n){
      yij = (x[j] + x[colp]) / 2;         /* x(i,j)+x(j,i) */
      y[j] = yij;
      y[colp] = yij;
    }
  }
}

/* ************************************************************
   PROCEDURE skewproj -- Y = (X-X')/2
   INPUT x, n - full n x n matrix x.
   OUTPUT y - on output, contains (x-x')/2
   ************************************************************ */
void skewproj(double *y, const double *x, const mwIndex n)
{
  mwIndex colp,i,j;
  double yij;

  /* ------------------------------------------------------------
     x points to x(:,i);     x+colp = x(:,j).
     ------------------------------------------------------------ */
  for(i = 0; i < n; x += n, y += n, i++){
    y[i] = 0.0;                         /* diagonal entry */
    for(colp = n + i, j=i+1; j<n; j++, colp += n){
      yij = (x[j] - x[colp]) / 2;         /* x(j,i) - x(i,j) */
      y[j] = yij;
      y[colp] = -yij;                   /* conjugate */
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
 mxArray *output_array[3], *Xk, *hXk;
 coneK cK;
 mwIndex k, nk, nkp1, nksqr, lendiag,lenud, lenfull, i,ii;
 double *lab,*q,*labk,*xwork;
 const double *x;

/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "psdeig requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "psdeig produces less output arguments");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Compute statistics based on cone K structure
   ------------------------------------------------------------ */
  lendiag = cK.rLen + cK.hLen;
  lenud = cK.rDim + cK.hDim;
  lenfull = cK.lpN + cK.qDim + lenud;
/* ------------------------------------------------------------
   Get input vector x
   ------------------------------------------------------------ */
  mxAssert(!mxIsSparse(X_IN), "x must be full (not sparse).");
  x = mxGetPr(X_IN);
  if(mxGetM(X_IN) * mxGetN(X_IN) != lenud){
    mxAssert(mxGetM(X_IN) * mxGetN(X_IN) == lenfull, "Size mismatch x");
    x += cK.lpN + cK.qDim;       /* point to PSD part */
  }

/* ------------------------------------------------------------
   Allocate output LAB(diag), eigvec Q(full for psd)
   ------------------------------------------------------------ */
  LAB_OUT = mxCreateDoubleMatrix(lendiag, (mwSize)1, mxREAL);
  lab = mxGetPr(LAB_OUT);
  if(nlhs > 1){
    Q_OUT = mxCreateDoubleMatrix(lenud, (mwSize)1, mxREAL);
    q = mxGetPr(Q_OUT);
  }
/* ------------------------------------------------------------
   Allocate working arrays:
   ------------------------------------------------------------ */
  Xk = mxCreateDoubleMatrix((mwSize)0,(mwSize)0,mxREAL);
  hXk = mxCreateDoubleMatrix((mwSize)0,(mwSize)0,mxCOMPLEX);
  xwork =(double *) mxCalloc(MAX(1,SQR(cK.rMaxn)+2*SQR(cK.hMaxn)),
                             sizeof(double));
/* ------------------------------------------------------------
   PSD: (I) LAB = eig(X)
   ------------------------------------------------------------ */
  if(nlhs < 2){
    for(k=0; k < cK.rsdpN; k++){                /* real symmetric */
      nk = cK.sdpNL[k];
      symproj(xwork,x,nk);              /* make it symmetric */
      mxSetM(Xk, nk);
      mxSetN(Xk, nk);
      mxSetPr(Xk, xwork);
      mexCallMATLAB(1, output_array, 1, &Xk, "eig");
      memcpy(lab, mxGetPr(output_array[0]), nk * sizeof(double));
/* ------------------------------------------------------------
   With mexCallMATLAB, we invoked the mexFunction "eig", which
   allocates a matrix struct *output_array[0], AND a block for the
   float data of that matrix.
   ==> mxDestroyArray() does not only free the float data, it
   also releases the matrix struct (and this is what we want).
   ------------------------------------------------------------ */
      mxDestroyArray(output_array[0]);
      lab += nk;  x += SQR(nk);
    }
/* ------------------------------------------------------------
   WARNING: Matlab's eig doesn't recognize Hermitian, hence VERY slow
   ------------------------------------------------------------ */
    for(; k < cK.sdpN; k++){                    /* complex Hermitian */
      nk = cK.sdpNL[k]; nksqr = SQR(nk);
      symproj(xwork,x,nk);              /* make it Hermitian */
      skewproj(xwork + nksqr,x+nksqr,nk);
      mxSetM(hXk, nk);
      mxSetN(hXk, nk);
      mxSetPr(hXk, xwork);
      mxSetPi(hXk, xwork + nksqr);     
      mexCallMATLAB(1, output_array, 1, &hXk, "eig");
      memcpy(lab, mxGetPr(output_array[0]), nk * sizeof(double));
      mxDestroyArray(output_array[0]);
      lab += nk;  x += 2 * nksqr;
    }
  } /* nlhs < 2 */
  else{
/* ------------------------------------------------------------
   SDP: (II) (Q,LAB) = eig(X)
   ------------------------------------------------------------ */
    for(k=0; k < cK.rsdpN; k++){                /* real symmetric */
      nk = cK.sdpNL[k];
      symproj(xwork,x,nk);                      /* make it symmetric */
      mxSetM(Xk, nk);
      mxSetN(Xk, nk);
      mxSetPr(Xk, xwork);
      mexCallMATLAB(2, output_array, 1, &Xk, "eig");
      nksqr = SQR(nk);                                  /* copy Q-matrix */
      memcpy(q, mxGetPr(output_array[0]), nksqr * sizeof(double));
      nkp1 = nk + 1;                                   /* copy diag(Lab) */
      labk = mxGetPr(output_array[1]);
      for(i = 0, ii = 0; i < nk; i++, ii += nkp1)
        lab[i] = labk[ii];
      mxDestroyArray(output_array[0]);
      mxDestroyArray(output_array[1]);
      lab += nk;  x += nksqr; q += nksqr;
    }
    for(; k < cK.sdpN; k++){                    /* complex Hermitian */
      nk = cK.sdpNL[k]; nksqr = SQR(nk);
      symproj(xwork,x,nk);                      /* make it Hermitian */
      skewproj(xwork + nksqr,x+nksqr,nk);
      mxSetM(hXk, nk);
      mxSetN(hXk, nk);
      mxSetPr(hXk, xwork);
      mxSetPi(hXk, xwork+nksqr);
#ifdef USE_SVD
      mexCallMATLAB(3, output_array, 1, &hXk, "svd");
#else
      mexCallMATLAB(2, output_array, 1, &hXk, "eig");
#endif
      memcpy(q, mxGetPr(output_array[0]), nksqr * sizeof(double));
      q += nksqr;
      if(mxIsComplex(output_array[0]))     /* if any imaginary part */
        memcpy(q, mxGetPi(output_array[0]), nksqr * sizeof(double));
      nkp1 = nk + 1;                              /* copy diag(Lab) */
      labk = mxGetPr(output_array[1]);
      for(i = 0, ii = 0; i < nk; i++, ii += nkp1)
        lab[i] = labk[ii];
      mxDestroyArray(output_array[0]);
      mxDestroyArray(output_array[1]);
#ifdef USE_SVD
      mxDestroyArray(output_array[2]);
#endif
      lab += nk;  x += 2 * nksqr; q += nksqr;
    }
  } /* [lab,q] = eigK */
/* ------------------------------------------------------------
   Release PSD-working arrays.
   ------------------------------------------------------------ */
  mxSetM(Xk,(mwSize)0); mxSetN(Xk,(mwSize)0); 
  mxSetPr(Xk, (double *) NULL);
  mxDestroyArray(Xk);
  mxSetM(hXk,(mwSize)0); mxSetN(hXk,(mwSize)0); 
  mxSetPr(hXk, (double *) NULL);   mxSetPi(hXk, (double *) NULL);
  mxDestroyArray(hXk);
  mxFree(xwork);
}
