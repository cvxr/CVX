/*
%                                      [zhi,zlo] = quadadd(xhi,xlo,y)
% QUADADD   Compute (zhi+zlo) = (xhi+xlo) + y.
%   x an z are in double-double format (quad precision)
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function [zhi,zlo] = quadadd(xhi,xlo,y)

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
#include "mex.h"
#include "blksdp.h"
#include <string.h>

#define ZHI_OUT myplhs[0]
#define ZLO_OUT myplhs[1]
#define NPAROUT 2

#define XHI_IN prhs[0]
#define XLO_IN prhs[1]
#define Y_IN prhs[2]
#define NPARIN 3

/* ************************************************************
   PROCEDURE quadadd
   RETURNS zhi.
   ************************************************************ */
double rquaddadd(double *pzlo, double xhi,double xlo, const double y)
{
  double zhi, zlo;
/* ------------------------------------------------------------
   |y| > |xhi| > |xlo|: Then fl(y+xlo) = y, thus add immediately to xhi.
   ------------------------------------------------------------ */
  if(fabs(y) > fabs(xhi)){
    zhi = y + xhi;
    xhi -= (zhi-y);         /* Residual xhi := xhi - fl((xhi+y)-y) */
    *pzlo = xlo + xhi;
    return zhi;
  }
/* ------------------------------------------------------------
   |xhi| > |y| > |xlo| or |xhi| > |xlo| > |y|
   ------------------------------------------------------------ */
  else{
    zlo = xlo + y;
    xlo -= (zlo - y);      /* If |xlo| > |y| then xlo := 0, otherwise resid. */
    zhi = (xhi + zlo);
    zlo -= (zhi - xhi);    /* resid. */
    *pzlo = xlo + zlo;   /* If |xlo| > |y| then zlo, which lost y-digits */
  }
  return zhi;
}

/* ============================================================
   MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
                 const int nrhs, const mxArray *prhs[])
{
  mxArray *myplhs[NPAROUT];
  mwSize m;
  mwIndex i;
  const double *xhi, *xlo, *y;
  double *zhi, *zlo;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "quadadd requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "quadadd produces less output arguments.");
/* ------------------------------------------------------------
   Get input vectors xhi,xlo,y.
   ------------------------------------------------------------ */
  m = mxGetM(XHI_IN);
  mxAssert(mxGetM(XLO_IN) == m, "xlo size mismatch.");
  mxAssert(mxGetM(Y_IN) == m, "y size mismatch.");
  xhi = mxGetPr(XHI_IN);
  xlo = mxGetPr(XLO_IN);
  y = mxGetPr(Y_IN);
/* ------------------------------------------------------------
   Create output vectors zhi(m), zlo(m).
   ------------------------------------------------------------ */
  ZHI_OUT = mxCreateDoubleMatrix(m,(mwSize)1,mxREAL);
  ZLO_OUT = mxCreateDoubleMatrix(m,(mwSize)1,mxREAL);
  zhi = mxGetPr(ZHI_OUT);
  zlo = mxGetPr(ZLO_OUT);
/* ------------------------------------------------------------
   Let (zhi,zlo) = quadadd(xhi,xlo, y)
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++)
    zhi[i] = rquaddadd(zlo+i, xhi[i],xlo[i],y[i]);
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  i = MAX(nlhs, 1);
  memcpy(plhs,myplhs, i * sizeof(mxArray *));
  for(; i < NPAROUT; i++)
    mxDestroyArray(myplhs[i]);
}
