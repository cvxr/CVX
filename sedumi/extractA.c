/* ************************************************************
%              Apart = extractA(At,Ajc,blk0,blk1,blkstart[,blkstart2])
% EXTRACTA  Fast alternative to
%  Apart = At(blkstart(1):blkstart(2)-1,:).
%  Instead of blkstart(2), it takes "blkstart2" (if supplied) or
%  size(At,1)+1 (if neither blkstart(2) nor blkstart2) are available.
%
%  Extract submatrix of
%  A with subscripts "blkstart(1):blkstart(2)-1" and has its nonzeros per column
%  bounded bij Ajc1:Ajc2. If Ajc1 = [] (Ajc2=[]) then start at column start
%  (end at column end) of A.
%
%  SEE ALSO partitA.
% ******************** INTERNAL FUNCTION OF SEDUMI ********************

function Apart = extractA(At,Ajc,blk0,blk1,blkstart[,blkstart2]) --

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

#define APART_OUT plhs[0]
#define NPAROUT 1

#define AT_IN prhs[0]
#define AJC_IN prhs[1]
#define BLK0_IN prhs[2]
#define BLK1_IN prhs[3]
#define BLKSTART_IN prhs[4]
#define NPARINMIN 5
#define BLKSTART2_IN prhs[5]
#define NPARIN 6

/* ************************************************************
   PROCEDURE extractA
   INPUT
     Ajc1, Ajc2, Air,Apr, m  - sparse N x m matrix, we only consider nonzeros
       in Air[Ajc1[k]:Ajc2[k]-1], k=0:m-1.
     ifirst - subscript of 1st entry, thus to be subtracted from A-subscripts
   OUTPUT
     Y - sparse n x m matrix, with sum(Ajc2-Ajc1) nonzeros.
   ************************************************************ */
void extractA(jcir Y, const mwIndex *Ajc1,const mwIndex *Ajc2,
              const mwIndex *Air, const double *Apr, const mwIndex ifirst,
              const mwIndex m)
{
  mwIndex i,j,inz,ajnnz;
  inz = 0;
  for(j = 0; j < m; j++){
    Y.jc[j] = inz;
    ajnnz = Ajc2[j]-Ajc1[j];
    memcpy(Y.pr + inz, Apr+Ajc1[j], ajnnz * sizeof(double));
    memcpy(Y.ir + inz, Air+Ajc1[j], ajnnz * sizeof(mwIndex));
    inz += ajnnz;
  }
/* ------------------------------------------------------------
   Close last column of Y, and subtract "ifirst" from its subscripts
   ------------------------------------------------------------ */
  Y.jc[m] = inz;
  for(i = 0; i < inz; i++)
    Y.ir[i] -= ifirst;
}

/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  jcir At, Apart;
  mwIndex i, m, njc, ifirst, n, ynnz, blk0,blk1;
  mwIndex *Ajc;
  const double *blkstartPr, *AjcPr;
  bool isblk0negative;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARINMIN, "extractA requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "extractA produces less output arguments.");
/* --------------------------------------------------
   GET inputs At, blkstart, Ajc1, Ajc2
   -------------------------------------------------- */
  mxAssert(mxIsSparse(AT_IN), "At must be a sparse matrix.");
  At.jc = mxGetJc(AT_IN);
  At.ir = mxGetIr(AT_IN);
  At.pr = mxGetPr(AT_IN);
  m = mxGetN(AT_IN);
  if(nrhs >=  NPARIN){
    n = (mwIndex) mxGetScalar(BLKSTART2_IN);
    ifirst = (mwIndex) mxGetScalar(BLKSTART_IN);
  }
  else if(mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN) ==2){
    blkstartPr = mxGetPr(BLKSTART_IN);
    ifirst = (mwIndex) blkstartPr[0];
    n = (mwIndex) blkstartPr[1];
  }
  else
    mxAssert(0==1, "blkstart size mismatch.");
  n -= (ifirst--);
  AjcPr = mxGetPr(AJC_IN);
  mxAssert(m == mxGetM(AJC_IN), "Ablkjc size mismatch.");
  njc = mxGetN(AJC_IN);
  blk0 = (mwIndex) mxGetScalar(BLK0_IN);           /* double to mwIndex */
  

  if(blk0>0){
    --blk0;                                /* Fortran to C */
    isblk0negative=0;}
  else
    isblk0negative=1;
  
 
  if(mxGetM(BLK1_IN) * mxGetN(BLK1_IN) != 1)
    blk1 = njc;                           /*default to end */
  else{
    blk1 = (mwIndex) mxGetScalar(BLK1_IN);   /* double to mwIndex (thus inf not allowed) */
    mxAssert(blk1>0,"");
    --blk1;                                /* Fortran to C */
  }
/* ------------------------------------------------------------
   Allocate mwIndex array blkstart Ajc(2*m)
   ------------------------------------------------------------ */
  Ajc = (mwIndex *) mxCalloc(MAX(2*m,1), sizeof(mwIndex));
/* ------------------------------------------------------------
   Convert Ajc from double to mwIndex:
   ------------------------------------------------------------ */
  mxAssert(blk0 < njc, "Ablkjc size mismatches blk0.");
  if(isblk0negative)
    memcpy(Ajc,At.jc,m*sizeof(mwIndex));          /* default: start of column */
  else
    for(i = 0; i < m; i++){                         /* to integers */
      Ajc[i] = (mwIndex) AjcPr[m*blk0 + i];
    }
  
  mxAssert(blk1 >= 0, "blk1 must be positive.");
  if(blk1 >= njc)
    memcpy(Ajc+m,At.jc+1,m*sizeof(mwIndex));      /* default: end of column */
  else
    for(i = 0; i < m; i++){                         /* to integers */
      Ajc[m+i] = (mwIndex) AjcPr[blk1*m + i];
    }
/* ------------------------------------------------------------
   Apart = sparse(n,m,ynnz);
   ------------------------------------------------------------ */
  ynnz = 0;
  for(i = 0; i < m; i++)
    ynnz += Ajc[m+i]-Ajc[i];              /* nnz(Apart(:,i)) */
  APART_OUT = mxCreateSparse(n,m, ynnz,mxREAL);
  Apart.jc = mxGetJc(APART_OUT);
  Apart.ir = mxGetIr(APART_OUT);
  Apart.pr = mxGetPr(APART_OUT);
/* ------------------------------------------------------------
   The real job:
   ------------------------------------------------------------ */
  extractA(Apart, Ajc,Ajc+m,At.ir,At.pr, ifirst, m);
/* ------------------------------------------------------------
   Release working array
   ------------------------------------------------------------ */
  mxFree(Ajc);
}
