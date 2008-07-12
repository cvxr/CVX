/*
%                                             y = givensrot(gjc,g,x,K)
% GIVENSROT
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = givensrot(gjc,g,x,K)

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
#include "givens.h"

#define Y_OUT plhs[0]
#define NPAROUT 1

#define GJC_IN prhs[0]
#define G_IN prhs[1]
#define X_IN prhs[2]
#define K_IN prhs[3]
#define NPARIN 4

/* ************************************************************
   PROCEDURE matgivens
   INPUT
     gjc, g - sequence of givens rotations
   UPDATED
     y - let Y := Q_g * Y.
   ************************************************************ */
void matgivens(double *y, const twodouble *g, const mwIndex *gjc, const mwIndex n)
{
  mwIndex j,k;

  for(j = 0; j < n; j++, y += n)       /* For all n columns of y */
/* ------------------------------------------------------------
   At step k, we apply m rotations involving rows k:k+m, m=gjc[k+1]-gjc[k].
   ------------------------------------------------------------ */
    for(k = 0; k < n-1; k++)
      givensrot(y+k, g+gjc[k], gjc[k+1]-gjc[k]);
}

/* complex case */
void prpimatgivens(double *y,double *ypi, const tridouble *g,
                   const mwIndex *gjc, const mwIndex n)
{
  mwIndex j,k;

  for(j = 0; j < n; j++, y += n, ypi += n)       /* For all n columns of y */
/* ------------------------------------------------------------
   At step k, we apply m rotations involving rows k:k+m, m=gjc[k+1]-gjc[k].
   ------------------------------------------------------------ */
    for(k = 0; k < n-1; k++)
      prpigivensrot(y+k,ypi+k, g+gjc[k], gjc[k+1]-gjc[k]);
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
  mwIndex inz, i, k, nk, nksqr, lenud, sdplen, gnnz;
  mwIndex *gjc, *iwork;
  const double *gjcPr;
  const double *g, *gk;
  double *y;
  coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "givensrot requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "givensrot generates less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  lenud = cK.rDim + cK.hDim;
  sdplen = cK.rLen + cK.hLen;
/* ------------------------------------------------------------
   Get inputs gjc,g,x
   ------------------------------------------------------------ */
  mxAssert(mxGetM(GJC_IN) * mxGetN(GJC_IN) == sdplen, "gjc size mismatch");
  gjcPr = mxGetPr(GJC_IN);
  g = (double *) mxGetPr(G_IN);
  gnnz = mxGetM(G_IN) * mxGetN(G_IN);
  mxAssert(mxGetM(X_IN) == lenud && mxGetN(X_IN) == 1, "x size mismatch");
/* ------------------------------------------------------------
   Allocate output y(lenud), and let y = x.
   ------------------------------------------------------------ */
  Y_OUT = mxCreateDoubleMatrix(lenud, (mwSize)1, mxREAL);
  y = mxGetPr(Y_OUT);
  memcpy(y, mxGetPr(X_IN), lenud * sizeof(double));
/* ------------------------------------------------------------
   Allocate working array iwork(sum(K.s))
   ------------------------------------------------------------ */
  iwork = (mwIndex *) mxCalloc(MAX(1,sdplen), sizeof(mwIndex));
/* ------------------------------------------------------------
   Convert gjcPr from float to mwIndex, and store in gjc:=iwork.
   ------------------------------------------------------------ */
  gjc = iwork;
  for(i = 0; i < sdplen; i++)
    gjc[i] = gjcPr[i];                 /* don't subtract 1: already C-style */
/* ------------------------------------------------------------
   The actual job is done here: U_NEW = Q(g) * U_OLD
   ------------------------------------------------------------ */
  inz = 0;
  for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
    nk = cK.sdpNL[k];
    nksqr = SQR(nk);
    mxAssert(inz + 2 * gjc[nk-1] <= gnnz, "g size mismatch");
    gk = g+inz;
    matgivens(y, (twodouble *) gk, gjc, nk);
    y += nksqr;
    inz += 2 * gjc[nk-1];        /* each rotation consists of 2 doubles */
    gjc += nk;
  }
  for(; k < cK.sdpN; k++){                       /* complex Hermitian */
    nk = cK.sdpNL[k];
    nksqr = SQR(nk);
    mxAssert(inz + 3 * gjc[nk-1] <= gnnz, "g size mismatch");
    gk = g+inz;
    prpimatgivens(y,y+nksqr, (tridouble *) gk, gjc, nk);
    nksqr += nksqr;
    y += nksqr;
    inz += 3 * gjc[nk-1];            /* each rotation consists of 3 doubles */
    gjc += nk;
  }
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(iwork);
}
