/* ************************************************************
%                                           perm = sortnnz(At,Ajc1,Ajc2)
% SORTNNZ  Sorts columns in At
%     in increasing order of nnzs; only the nnzs between Ajc1 and Ajc2
%     are considered for each column. If Ajc1 or Ajc2 is empty, we use
%     the start or end of the columns in At.
%
%  SEE ALSO partitA.
% ******************** INTERNAL FUNCTION OF SEDUMI ********************

function perm = sortnnz(At,Ajc1,Ajc2)

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

   ************************************************************ */
#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define PERM_OUT plhs[0]
#define NPAROUT 1

#define AT_IN prhs[0]
#define AJC1_IN prhs[1]
#define AJC2_IN prhs[2]
#define NPARIN 3

/* ************************************************************
   PROCEDURE sortnnz - sort "jc2-jc1" in increasing order
   INPUT
     jc1, jc2 - length m arrays
     m  - order
   OUTPUT
     y - length m keyint-array.
   ************************************************************ */
void sortnnz(keyint *y, const mwIndex *jc1, const mwIndex *jc2, const mwSize m)
{
  keyint ki;
  for(ki.k = 0; ki.k < m; ki.k++){
    ki.i = jc2[ki.k] - jc1[ki.k];
    y[ki.k] = ki;
  }
  kiqsort(y,m);
}


/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  mwSize i, m, njc1,njc2;
  mwIndex *Ajc;
  const double *Ajc1Pr, *Ajc2Pr;
  double *permPr;
  keyint *y;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "sortnnz requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "sortnnz produces less output arguments.");
/* --------------------------------------------------
   GET inputs At, Ajc1, Ajc2
   -------------------------------------------------- */
  mxAssert(mxIsSparse(AT_IN), "At must be a sparse matrix.");
  m = mxGetN(AT_IN);
  Ajc1Pr = mxGetPr(AJC1_IN);
  njc1 = mxGetM(AJC1_IN) * mxGetN(AJC1_IN);
  Ajc2Pr = mxGetPr(AJC2_IN);
  njc2 = mxGetM(AJC2_IN) * mxGetN(AJC2_IN);
/* ------------------------------------------------------------
   Allocate int array Ajc(2*m)
   keyint array y(m)
   ------------------------------------------------------------ */
  Ajc = (mwIndex *) mxCalloc(MAX(2*m,1), sizeof(mwIndex));
  y = (keyint *) mxCalloc(MAX(m,1), sizeof(keyint));
/* ------------------------------------------------------------
   Convert Ajc from double to int:
   ------------------------------------------------------------ */
  if(njc1 == 0)
    memcpy(Ajc,mxGetJc(AT_IN),m*sizeof(mwIndex));  /* default: start of column */
  else {
    mxAssert(njc1 >= m, "Ajc1 size mismatch.");
    for(i = 0; i < m; i++){                         /* to integers */
      Ajc[i] = Ajc1Pr[i];
    }
  }
  if(njc2 == 0)
    memcpy(Ajc+m,mxGetJc(AT_IN)+1,m*sizeof(mwIndex));  /* default: end of column */
  else {
    mxAssert(njc2 >= m, "Ajc2 size mismatch.");
    for(i = 0; i < m; i++){                         /* to integers */
      Ajc[m+i] = Ajc2Pr[i];
    }
  }
/* ------------------------------------------------------------
   Create output perm(m)
   ------------------------------------------------------------ */
  PERM_OUT = mxCreateDoubleMatrix(m, (mwSize)1, mxREAL);
  permPr = mxGetPr(PERM_OUT);
/* ------------------------------------------------------------
   The real job:
   ------------------------------------------------------------ */
  sortnnz(y, Ajc, Ajc+m, m);
  for(i = 0; i < m; i++)
    permPr[i] = (y+i)->k + 1.0;
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(y);
  mxFree(Ajc);
}
