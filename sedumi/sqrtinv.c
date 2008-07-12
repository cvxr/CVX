/*
%                                                  y = sqrtinv(q,vlab,K)
% SQRTINV  Computes for PSD-cone, y = (Q / diag(sqrt(vlab)))', so that
%   Y'*Y = inv(Q * diag(vlab) * Q').
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = sqrtinv(q,vlab,K)

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

#define Q_IN prhs[0]
#define V_IN prhs[1]
#define K_IN prhs[2]

#define NPARIN 3

/* ************************************************************
   PROCEDURE qdivv - Computes Y = (Q / diag(sqrt(vlab)))'
   ************************************************************ */
void qdivv(double *y, const double *q,const double *v,const mwIndex n)
{
  mwIndex i,j,jcol,inz;
  double svi;

  for(inz = 0, i = 0; i < n; i++){
    svi = sqrt(v[i]);
    for(j = 0, jcol = i; j < n; j++, jcol += n)
      y[jcol] = q[inz++] / svi;              /* y(i,j) = q(j,i) / svi */
  }
}

/* complex case */
void prpiqdivv(double *y,double *ypi, const double *q,const double *qpi,
               const double *v,const mwIndex n)
{
  mwIndex i,j,jcol,inz;
  double svi;

  for(inz = 0, i = 0; i < n; i++){
    svi = sqrt(v[i]);
    for(j = 0, jcol = i; j < n; j++, jcol += n){
/* y(i,j) = conj(q(j,i)) / svi */
      y[jcol] = q[inz] / svi;
      ypi[jcol] = -qpi[inz] / svi;                    /* conjugate */
      inz++;
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
  mwIndex i,k, nk, nksqr, lenud, lendiag, diagskip;
  double *y;
  const double *q,*v;
  coneK cK;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "sqrtinv requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "sqrtinv generates less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  lenud = cK.rDim + cK.hDim;
  diagskip = cK.lpN + 2 * cK.lorN;         /* diag for LP and Lorentz */
  lendiag = diagskip + cK.rLen + cK.hLen;
/* ------------------------------------------------------------
   Get inputs v, q
   ------------------------------------------------------------ */
  mxAssert(mxGetM(V_IN) * mxGetN(V_IN) == lendiag, "v size mismatch");
  mxAssert(mxGetM(Q_IN) * mxGetN(Q_IN) == lenud, "q size mismatch");
  q = mxGetPr(Q_IN);
  v = mxGetPr(V_IN) + diagskip;    /* skip LP*/
/* ------------------------------------------------------------
   Allocate output y
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(lenud, (mwSize)1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   The actual job is done here:
   ------------------------------------------------------------ */
  for(k = 0; k < cK.rsdpN; k++){                /* PSD: real symmetric */
    nk = cK.sdpNL[k];
    qdivv(y, q,v,nk);
    nksqr = SQR(nk);
    y += nksqr; q += nksqr; v += nk;
  }
  for(; k < cK.sdpN; k++){                    /* complex Hermitian */
    nk = cK.sdpNL[k];
    nksqr = SQR(nk);
    prpiqdivv(y,y+nksqr, q,q+nksqr, v, nk);
    nksqr += nksqr;
    y += nksqr; q += nksqr; v += nk;
  }
}
