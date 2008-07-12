/*
%                                                y = qblkmul(mu,d,blkstart)
% QBLKMUL  yields length(y)=blkstart(end)-blkstart(1) vector with
%    y[k] = mu(k) * d[k]; the blocks d[k] are partitioned by blkstart.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = qblkmul(mu,d,blkstart)

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
#define BLKSTART_IN prhs[2]
#define NPARIN 3

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  mwIndex i,j,nblk,k,nk,qDim;
  double *y;
  const double *d, *mu, *blkstartPr;
  mwIndex *blkstart;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "qblkmul requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "qblkmul generates 1 output argument.");
/* ------------------------------------------------------------
   Get inputs d, mu, blkstart
   ------------------------------------------------------------ */
  d = mxGetPr(D_IN);
  qDim = mxGetM(D_IN) * mxGetN(D_IN);
  nblk = mxGetM(MU_IN) * mxGetN(MU_IN);
  mu = mxGetPr(MU_IN);
  blkstartPr = mxGetPr(BLKSTART_IN);
  mxAssert(nblk == mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN) - 1, "blkstart size mismatch.");
/* ------------------------------------------------------------
   Allocate mwIndex working array blkstart(nblk+1).
   ------------------------------------------------------------ */
  blkstart = (mwIndex *) mxCalloc(nblk + 1, sizeof(mwIndex));
/* ------------------------------------------------------------
   Convert Fortran double to C mwIndex
   ------------------------------------------------------------ */
  for(i = 0; i <= nblk; i++){
    j = blkstartPr[i];             /* double to mwIndex */
    blkstart[i] = --j;
  }
/* ------------------------------------------------------------
   Let d point to Lorentz norm-bound
   ------------------------------------------------------------ */
  if(qDim != blkstart[nblk] - blkstart[0]){
    if(qDim == nblk + blkstart[nblk] - blkstart[0]){
      d += nblk;                   /* Point to Lorentz norm-bound */
    }
    else {
      mxAssert(qDim >= blkstart[nblk], "d size mismatch.");
      d += blkstart[0];                   /* Point to Lorentz norm-bound */
    }
  }
  qDim = blkstart[nblk] - blkstart[0];
/* ------------------------------------------------------------
   Allocate output y(qDim)
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(qDim, (mwSize)1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   LORENTZ: yk = mu(k) * d[k]
   ------------------------------------------------------------ */
  for(k = 0; k < nblk; k++){
    nk = blkstart[k+1] - blkstart[k];
    scalarmul(y,mu[k],d,nk);
    y += nk;
    d += nk;
  }
}
