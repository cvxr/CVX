/*
 smult = adenscale(dense, d, qblkstart)

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

#define SMULT_OUT plhs[0]
#define NPAROUT 1

#define DENSE_IN prhs[0]
#define D_IN prhs[1]
#define BLKSTART_IN prhs[2]
#define NPARIN 3

/* ************************************************************
   PROCEDURE adenscale - Computed Lorentz norm-bound part of smult such that
     AP(d)A' = ADA + Ad*diag(smult)*Ad'. Thus, smult(j) = det(k) if j
     belongs to Lorentz block k.
   INPUT
     detd     - length |K.q| vector of det(dk)_k for Lorentz; use its sqrt.
     dencols  - length nden array
     q        - length nq array
     blkend   - length nq array, with 1-beyond subscript for each Lorentz
        block listed in q.
     nq       - number of (removed) dense Lorentz blocks
     nden     - number of dense columns, nden >= nl + nq.
   OUTPUT
     smult - length nden vector
   ************************************************************ */
void adenscale(double *smult, const double *detd, const mwIndex *dencols,
               const mwIndex *q, const mwIndex *blkend, const mwIndex nq, const mwIndex nden)
{
  mwIndex j,k;
  double detdk;
/* ------------------------------------------------------------
   LORENTZ norm-bound, detd(q(k)) while dencols(j)<blkend(k)
   ------------------------------------------------------------ */
  j= 0;
  for(k = 0; k < nq; k++){
    detdk = detd[q[k]];
    while(j < nden){
      if(dencols[j] >= blkend[k])
        break;
      smult[j++] = detdk;
    }
  }
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction( int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  const mxArray *MY_FIELD;
  mwIndex i, j, nden, nl, nq, lorN;
  mwIndex *q, *dencols, *blkend;
  const double *qPr, *dencolsPr, *detd, *blkstartPr;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "adenscale requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "adenscale produces less output arguments");
/* ------------------------------------------------------------
   DISASSEMBLE dense structure: dense.{l,cols,q}
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(DENSE_IN), "dense should be a structure.");
  MY_FIELD = mxGetField(DENSE_IN,(mwIndex)0,"l");        /* dense.l */ 
  mxAssert( MY_FIELD != NULL, "Missing field dense.l.");
  nl = (mwIndex) mxGetScalar(MY_FIELD);                           /* double to mwIndex */
  MY_FIELD = mxGetField(DENSE_IN,(mwIndex)0,"q");          /* dense.q */
  mxAssert( MY_FIELD != NULL, "Missing field dense.q.");
  nq = (mwIndex) (mxGetM(MY_FIELD) * mxGetN(MY_FIELD));
  qPr = mxGetPr(MY_FIELD);
  MY_FIELD = mxGetField(DENSE_IN,(mwIndex)0,"cols");            /* dense.cols */
  mxAssert( MY_FIELD != NULL, "Missing field dense.cols.");
  nden = (mwIndex) (mxGetM(MY_FIELD) * mxGetN(MY_FIELD) - (nq+nl));  /* jump to q-norm */
  mxAssert(nden >= 0, "dense.cols size mismatch.");
  dencolsPr = mxGetPr(MY_FIELD) + (nq+nl);               /* jump to q-norm */
/* ------------------------------------------------------------
   Disassemble structure d.{det}
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(D_IN), "d should be a structure.");
  MY_FIELD = mxGetField(D_IN,(mwIndex)0,"det");           /* d.det */
  mxAssert( MY_FIELD != NULL, "Missing field d.det.");
  detd = mxGetPr(MY_FIELD);
  lorN = (mwIndex) (mxGetM(MY_FIELD) * mxGetN(MY_FIELD));
/* ------------------------------------------------------------
   Get INPUTS blkstart
   ------------------------------------------------------------ */
  blkstartPr = mxGetPr(BLKSTART_IN);
  mxAssert(lorN +1 == mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN), "blkstart size mismatch");
/* ------------------------------------------------------------
   Create working arrays q(nq), blkend(nq), dencols(nden - (nl+nq))
   ------------------------------------------------------------ */
  q       = (mwIndex *) mxCalloc(MAX(nq,1), sizeof(mwIndex));
  blkend  = (mwIndex *) mxCalloc(MAX(nq,1), sizeof(mwIndex));
  dencols = (mwIndex *) mxCalloc(MAX(1,nden), sizeof(mwIndex));
/* ------------------------------------------------------------
   Convert to integer C-style: q, blkend, dencols.
   ------------------------------------------------------------ */
  for(i = 0; i < nq; i++){
    j = (mwIndex) qPr[i];
    mxAssert(j>0,"");
    q[i] = --j;
  }
  for(i = 0; i < nq; i++){
    j = (mwIndex) blkstartPr[q[i] + 1];        /* F-double to C-mwIndex */
    mxAssert(j>0,"");    
    blkend[i] = --j;
  }
  for(i = 0; i < nden; i++){
    j = (mwIndex) dencolsPr[i];
    mxAssert(j>0,"");    
    dencols[i] = --j;
  }
/* ------------------------------------------------------------
   Create output: smult(nden,1)
   ------------------------------------------------------------ */
  SMULT_OUT = mxCreateDoubleMatrix(nden, (mwSize)1, mxREAL);
  adenscale(mxGetPr(SMULT_OUT),detd, dencols,q, blkend, nq,nden);
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(dencols);
  mxFree(blkend);
  mxFree(q);
}
