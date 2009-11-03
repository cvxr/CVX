/*
    y = mJdetd(detd,K)
    yields length(y)=sum(K.q) vector with y[k] = -detd(k)*J,
    where J = [1, 0'; 0, -I].

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
#define DETD_IN prhs[0]
#define K_IN prhs[1]
#define NPARIN 2

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   y = mJdetd(detd,K)
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
 mwIndex i,nexti,k;
 double detdk;
 double *y;
 const double *detd;
 coneK cK;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
 if(nrhs < NPARIN)
   mexErrMsgTxt("mJdetd requires more input arguments.");
 if (nlhs > NPAROUT)
   mexErrMsgTxt("mJdetd generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
 conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get input detd
   ------------------------------------------------------------ */
 if(mxGetM(DETD_IN) * mxGetN(DETD_IN) != cK.lorN)
   mexErrMsgTxt("detd size mismatch");
 detd = mxGetPr(DETD_IN);
/* ------------------------------------------------------------
   Allocate output y(qDim)
   ------------------------------------------------------------ */
 Y_OUT =  mxCreateDoubleMatrix(cK.qDim, 1, mxREAL);
 y = mxGetPr(Y_OUT);
/* ------------------------------------------------------------
   LORENTZ: yk = [-detd(k); detd(k)* ones(nk-1,1)]
   ------------------------------------------------------------ */
 i = 0;
 nexti = 0;
 for(k = 0; k < cK.lorN; k++){
   nexti += cK.lorNL[k];
   detdk = detd[k];
   y[i] = -detdk;
   for(++i; i < nexti; i++)
     y[i] = detdk;
 }
}
