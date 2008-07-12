/*
   x = psdframeit(lab,frms,K);

   Computes x = FRM*lab.
   The frame "frms" is a product-form Householder reflection. 

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
#include "reflect.h"

#define X_OUT plhs[0]
#define NPAROUT 1

#define LAB_IN prhs[0]
#define FRMS_IN prhs[1]
#define K_IN prhs[2]
#define NPARIN 3

/* ************************************************************
   PROCEDURE psdframeit
   INPUT
     frms - lenud+hlen eigenvector basis/ orthogonal matrix Q,
       in Householder form.
     lab - sdplen vector of eigenvalues
     sdpNL sdpN-array: order of each PSD block, i.e. K.s.
     rsdpN number of real symmetric blocks: 0<= rsdpN <= sdpN.
     sdpN number of PSD blocks, i.e. length(K.s).
   UPDATED
     x - lenud vector. On input zeros. On output, x = Q*LAB*Q'.
   WORKING ARRAY:
     fwork length max(rmaxn,2*hmaxn).
   ************************************************************ */
void psdframeit(double *x, const double *frms, const double *lab,
                const mwIndex *sdpNL,const mwIndex rsdpN,const mwIndex sdpN,
                double *fwork)
{
  mwIndex k,nk,nksqr, i,inz;
  const double *beta;
  for(k = 0; k < rsdpN; k++){                /* real symmetric */
    nk = sdpNL[k];
    nksqr = SQR(nk);
    for(i = 0, inz = 0; i < nk; i++, inz += nk+1)  /* Let X = diag(lab) */
      x[inz] = lab[i];
    qtxq(x, frms + nksqr - nk, frms, nk, fwork);
    tril2sym(x,nk);
    x += nksqr; frms += nksqr;
    lab += nk;
  }
/* ------------------------------------------------------------
   For complex Hermitian, we have X = Qb' * diag(lab) * Qb,
   where Qb = Q_1 * ... * Q_{n-1} * diag(q(:,n)),
   with q(:,n) a complex sign-vector.
   Format: [RE q(:,1:n-1), RE q(:,n), IM q(:,1:n-1), IM q(:,n), beta].
   ------------------------------------------------------------ */
  for(; k < sdpN; k++){                    /* complex Hermitian */
    nk = sdpNL[k];
    nksqr = SQR(nk);
    for(i = 0, inz = 0; i < nk; i++, inz += nk+1)  /* Let X = diag(lab) */
      x[inz] = lab[i];
    beta = frms + 2*nksqr;                   /* beta = frms(:,2*n+1) */
    prpiqtxq(x,x+nksqr, beta, frms,frms+nksqr, nk, fwork);
    tril2herm(x,x+nksqr, nk);
    nksqr += nksqr;
    x += nksqr; frms += nksqr + nk;          /* skip also beta. */
    lab += nk;
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
  mwIndex i,lendiag, lenfull, lenud,qsize;
  double *x, *fwork;
  const double *lab,*frms;
  mwIndex *sdpNL;
  coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "psdframeit requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "psdframeit generates less output arguments.");
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
   Get inputs lab,frms
   ------------------------------------------------------------ */
  lab = mxGetPr(LAB_IN);
  if(mxGetM(LAB_IN) * mxGetN(LAB_IN) != cK.rLen + cK.hLen){
    mxAssert(mxGetM(LAB_IN) * mxGetN(LAB_IN) == lendiag, "lab size mismatch");
    lab += cK.lpN + 2 * cK.lorN;
  }
  mxAssert(mxGetM(FRMS_IN) * mxGetN(FRMS_IN) == qsize, "frms size mismatch");
  frms = mxGetPr(FRMS_IN);
/* ------------------------------------------------------------
   Allocate output x
   ------------------------------------------------------------ */
  X_OUT =  mxCreateDoubleMatrix(lenud, (mwSize)1, mxREAL);
  x = mxGetPr(X_OUT);
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
    sdpNL[i] = cK.sdpNL[i];
/* ------------------------------------------------------------
   PSD: X = Qb' * diag(lab) * Qb,          Qb = Q_1*Q_2*..*Q_{n-1}
   where Q_k = I-ck*ck'/betak is an elementary Householder reflection.
   Format: frms = [c1, .., c_{n-1}, beta].
   VERY INEFFICIENT !!!!
   ------------------------------------------------------------ */
  psdframeit(x, frms,lab,sdpNL,cK.rsdpN,cK.sdpN,fwork);
/* ------------------------------------------------------------
   Release working array
   ------------------------------------------------------------ */
  mxFree(sdpNL);
  mxFree(fwork);
}
