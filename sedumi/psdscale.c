/*
%                                           y = psdscale(ud,x,K [,transp])
% PSDSCALE  Computes length lenud (=sum(K.s.^2)) vector y.
%   !transp (default) then y[k] = vec(Ldk' * Xk * Ldk)
%   transp == 1 then y[k] = vec(Udk' * Xk * Udk)
%   Uses pivot ordering ud.perm if available and nonempty.
%
%   SEE ALSO scaleK, factorK.
% **********  INTERNAL FUNCTION OF SEDUMI **********
function y = psdscale(ud,x,K [,transp])

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

/*    y = psdscale(ud,x,K [,transp]) */

#define Y_OUT plhs[0]
#define NPAROUT 1

#define UD_IN prhs[0]
#define X_IN prhs[1]
#define K_IN prhs[2]
#define NPARINMIN 3
#define TRANSP_IN prhs[3]
#define NPARIN 4

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  const mxArray *UD_FIELD;
  mwIndex lenfull, lenud, sdplen, fwsiz, i,k;
  double *fwork, *y,*permPr;
  const double *x,*ud;
  mwIndex *perm, *iwork;
  coneK cK;
  bool use_pivot, transp;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  if(nrhs < NPARIN){
    transp = 0;
    mxAssert(nrhs >= NPARINMIN, "psdscale requires more input arguments.");
  }
  else
    transp = (bool)mxGetScalar(TRANSP_IN);
  mxAssert(nlhs <= NPAROUT, "psdscale generates less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  lenud = cK.rDim + cK.hDim;
  lenfull = cK.lpN +  cK.qDim + lenud;
  sdplen = cK.rLen + cK.hLen;
/* ------------------------------------------------------------
   Get scale data: ud.{u,perm}.
   ------------------------------------------------------------ */
  if(!mxIsStruct(UD_IN)){
    mxAssert(mxGetM(UD_IN) * mxGetN(UD_IN) == lenud, "ud size mismatch."); /* ud is vector */
    ud = mxGetPr(UD_IN);
    use_pivot = 0;
  }
  else{                         /* ud is structure */
    UD_FIELD = mxGetField(UD_IN,(mwIndex)0,"u");
      mxAssert( UD_FIELD!= NULL,  "Field ud.u missing.");    /* ud.u */
    mxAssert(mxGetM(UD_FIELD) * mxGetN(UD_FIELD) == lenud, "ud.u size mismatch.");
    ud = mxGetPr(UD_FIELD);
    UD_FIELD = mxGetField(UD_IN,(mwIndex)0,"perm");                       /* ud.perm */
    if((use_pivot = (UD_FIELD != NULL))){
      if(mxGetM(UD_FIELD) * mxGetN(UD_FIELD) == sdplen)
        permPr = mxGetPr(UD_FIELD);
      else {
        mxAssert(mxGetM(UD_FIELD) * mxGetN(UD_FIELD) == 0, "ud.perm size mismatch");
        use_pivot = 0;
      }
    }
  }
/* ------------------------------------------------------------
   Get input x
   ------------------------------------------------------------ */
  mxAssert(!mxIsSparse(X_IN), "Sparse x not supported by this version of psdscale.");
  x = mxGetPr(X_IN);
/* ------------------------------------------------------------
   Validate x-input, and let x point to PSD part.
   ------------------------------------------------------------ */
  if(mxGetM(X_IN) * mxGetN(X_IN) == lenfull)
    x += cK.lpN + cK.qDim;
  else mxAssert(mxGetM(X_IN) * mxGetN(X_IN) == lenud, "size x mismatch.");
/* ------------------------------------------------------------
   Allocate output Y(lenud)
   ------------------------------------------------------------ */
  Y_OUT =  mxCreateDoubleMatrix(lenud, (mwSize)1, mxREAL);
  y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   Allocate fwork 2 * [ max(cK.rMaxn^2, 2*cK.hMaxn^2) ]
   iwork = mwIndex(sdplen)
   ------------------------------------------------------------ */
  fwsiz = MAX(SQR(cK.rMaxn),2*SQR(cK.hMaxn));
  fwork = (double *) mxCalloc( MAX(1,2 * fwsiz), sizeof(double));
  iwork = (mwIndex *) mxCalloc( MAX(1, sdplen), sizeof(mwIndex) );
/* ------------------------------------------------------------
   Convert Fortran to C-style in perm:
   ------------------------------------------------------------ */
  if(use_pivot){
    perm = iwork;
    for(k = 0; k < sdplen; k++){
      i = permPr[k];
      perm[k] = --i;
    }
  }
  else
    perm = (mwIndex *) NULL;
/* ------------------------------------------------------------
   The actual job is done here:.
   ------------------------------------------------------------ */
  psdscaleK(y, ud, perm, x, cK, transp, fwork);
/* ------------------------------------------------------------
   Release working arrays.
   ------------------------------------------------------------ */
  mxFree(fwork);
  mxFree(iwork);
}
