/*
%                                              y = qreshape(x,flag, K)
% QRESHAPE  Reshuffles entries associated with Lorentz blocks.
%   If flag = 0 then y = [x1 for each block; x2 for each block]
%   If flag = 1 then y = [x block 1; x block 2; etc], etc
%   Thus, x = qreshape(qreshape(x,0,K),1,K).
%
% SEE ALSO sedumi
% ******************** INTERNAL FUNCTION OF SEDUMI ********************

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

/* y=qreshape(x,flag,K)
  flag = 0 -> y = [x1 for each block; x2 for each block]
  flag = 1 -> y = [x block 1; x block 2; etc], etc
 thus, x = qreshape(qreshape(x,0,K),1,K).
*/

#include <string.h>
#include <math.h>      /* floor and log */
#include "mex.h"
#include "blksdp.h"

#define Y_OUT plhs[0]
#define NPAROUT 1
#define X_IN prhs[0]
#define FLAG_IN prhs[1]
#define K_IN prhs[2]
#define NPARIN 3

void intadd(mwIndex *x, const mwIndex y, const mwIndex n)
{
  mwIndex i;
  for(i = 0; i < n; i++)
    x[i] += y;
}

/* ************************************************************
   PROCEDURE qreshape0 - Transforms from block-wise to [x1(lorN); x2-blocks].
   INPUT
     blks - lorN array; blks[k]:blk[k+1]-1 are the subscripts of Lorentz
        block k=0:lorN-2.
     lorN - number of Lorentz blocks
   UPDATED
     y - entries in y(blk[0]:blks[lorN-1]) are affected (reshuffled)
   WORK
     fwork - length lorN-1 vector.
   ************************************************************ */
void qreshape0(double *y, const mwIndex *blks, const mwIndex lorN, double *fwork)
{
  mwIndex i,j,k;
  if(lorN <=1)
    return;             /* Nothing to do if only 1 block */
/* ------------------------------------------------------------
   Save y(2:lorN).
   Then make indices blks[0]+1:blks[0]+lorN-1 valid into fwork.
   ------------------------------------------------------------ */
  mxAssert(lorN >= 2,"");
  memcpy(fwork, y+blks[0] + 1, (lorN-1) * sizeof(double));
  fwork -= blks[0] + 1;    /* Make blks[0]+1 the 1st index */
/* ------------------------------------------------------------
   Let y(2:lorN) = y(iTr(2:end)), where iTr = cumsum([1 lorNL]).
   ------------------------------------------------------------ */
  for(k = 1; k < lorN; k++)
    y[blks[0] + k] = y[blks[k]];
/* ------------------------------------------------------------
   For each block k = 1:lorN-1, shift y2[lorN-k] k positions downwards,
   until we hit index "j". The data x(blks[0]+1:j-1) is in fwork.
   ------------------------------------------------------------ */
  j = blks[0] + lorN; /* 1st index beyond those saved in fwork*/
  for(k = 1; k < lorN; k++){
    if((i = blks[lorN-k-1] + 1) < j)
      break;
    if(i < blks[lorN-k])                   /* should be, since min(K.q)>=3 */
      memmove(y+i+k,y+i, (blks[lorN-k]-i) * sizeof(double));
  }
/* ------------------------------------------------------------
   Move block lorN-k, with part in fwork, part in y
   ------------------------------------------------------------ */
  if(k < lorN){
    if(blks[lorN-k] > j){        /* Still part in y ? */
      memmove(y+j+k,y+j, (blks[lorN-k]-j) * sizeof(double));
      mxAssert(i < j,"");
      memcpy(y+i+k,fwork+i, (j-i) * sizeof(double));
      k++;   /* finished this block */
    }
/* ------------------------------------------------------------
   Copy remaining blocks lorN-k, k < lorN, which are completely in fwork.
   ------------------------------------------------------------ */
    for(; k < lorN; k++){
      i = blks[lorN-k-1] + 1;
      if(i < blks[lorN-k])                   /* should be, since min(K.q)>=3 */
        memcpy(y+i+k,fwork+i, (blks[lorN-k]-i) * sizeof(double));
    }
  }
}

/* ************************************************************
   PROCEDURE spqreshape0 - Transforms from block-wise to [x1(lorN); x2-blocks].
   INPUT
     ynnz - nnz(y)
     blks - lorN+1 array; blks[k]:blk[k+1]-1 are the subscripts of Lorentz
        block k = 0:lorN-2.
     lorN - number of Lorentz blocks
     iwsize - length of iwork, should be 2*(1+lorN).
   UPDATED
     yir, ypr - entries with subscripts in y(blk[0]:blks[lorN]) are
       affected (reshuffled)
   WORK
     cfound - length lorN character working array.
     iwork - length iwsize = 2*(1+lorN).
     fwork - length lorN vector.
   ************************************************************ */
void spqreshape0(mwIndex *yir, double *ypr, mwIndex ynnz, const mwIndex *blks,
                 const mwIndex lorN, mwIndex iwsize, bool *cfound,
                 mwIndex *iwork, double *fwork)
{
  mwIndex inz,k,knz, blknnz;
  mwIndex *ipos;

  if(lorN <=1)
    return;             /* Nothing to do if only 1 block */
/* ------------------------------------------------------------
   Partition WORKING ARRAY:
   ipos(2+lorN), iwork.
   ------------------------------------------------------------ */
  ipos = iwork;
  iwork += 2 + lorN;
  iwsize -= 2 + lorN;
/* ------------------------------------------------------------
   Search for 1st and last Lorentz nonzero, let (yir,ypr) point
   to 1st lorentz nonzero, and let ynnz denote the number of them.
   ------------------------------------------------------------ */
  inz = 0;
  intbsearch(&inz,yir,ynnz,blks[0]);
  yir += inz;
  ypr += inz;
  ynnz -= inz;
  inz = 0;
  intbsearch(&inz,yir,ynnz,blks[lorN]);
  ynnz = inz;
  if(lorN <= ynnz){
/* ------------------------------------------------------------
   CASE 1: more nonzeros than blocks:
   Locate blks(0:lorN) in yir(1:xnnz)
   We'll have yir[ipos[k+1]-1] < blks[k] <= yir[ipos[k+1]],
   SO: yir[ipos[k+1]] is 1st nz belonging to block k OR LATER.
   NB: iwork(blknnz) lists nz-blocks, ipos(blknnz) points to 1st nonzero
     of that nonzero block in yir.
   ------------------------------------------------------------ */
    if(intmbsearch(ipos, cfound, yir, ynnz, blks, lorN, iwork, iwsize) != 0)
      mexErrMsgTxt("Out of working space");
    knz = 0;
    for(k = 1; k <= lorN; k++)
      if(ipos[k] < ipos[k+1]){
        iwork[knz] = k-1;           /* Store in C-style: 0:lorN-1 */
        ipos[knz] = ipos[k];        /* start of block k in yir */
        ++knz;
      }
    ipos[knz] = ipos[lorN+1];
    blknnz = knz;
    for(knz = 0; knz < blknnz; knz++)
      cfound[knz] = cfound[iwork[knz]];
  }
  else{
/* ------------------------------------------------------------
   CASE 2: more blocks than nonzeros:
   Locate nonzeros in list of blkstarts.  If cfound[inz] then yir[inz]
   is 1st of block ipos[inz+1], otherwise it is in block ipos[inz+1]-1.
     We'll have blks[ipos[inz+1]-1] < yir[inz] <= blks[ipos[inz+1]],
   so yir[inz] belongs to block ipos[inz+1]-1 if !cfound,
   and is tr-elt of ipos[inz+1] if cfound.
   NB: iwork(blknnz) lists nz-blocks, ipos(blknnz) points to 1st nonzero
     of that nonzero block in yir.
   ------------------------------------------------------------ */
    if(intmbsearch(ipos, cfound, blks, lorN, yir, ynnz, iwork, iwsize) != 0)
      mexErrMsgTxt("Out of working space");
    knz = 0;
    for(inz = 0; inz < ynnz; inz++) {
      if(cfound[inz]){                 /* block opens on trace */
        k = ipos[inz+1];
        iwork[knz] = k;
        ipos[knz] = inz;               /* start of block k in yir */
        cfound[knz++] = 1;
      }
      else if(knz==0||ipos[inz+1]>k+1){      /* block opens somewhere later */
        k = ipos[inz+1] - 1;
        iwork[knz] = k;
        ipos[knz] = inz;               /* start of block k in yir */
        cfound[knz++] = 0;
      }
    }
    ipos[knz] = ynnz;
    blknnz = knz;
  }
/* ------------------------------------------------------------
   Update all subscripts of norm-bound entries: += lorN-(k+1).
   ------------------------------------------------------------ */
  for(knz = 0; knz < blknnz; knz++){
    k = iwork[knz];
    if(++k < lorN){
      inz = ipos[knz] + cfound[knz];
      intadd(yir + inz, lorN-k, ipos[knz+1] - inz);
    }
  }
/* ------------------------------------------------------------
   Condense to a list of nonzero tr-entries
   ------------------------------------------------------------ */
  for(knz = 0; knz < blknnz; knz++)
    if(!cfound[knz])
      break;
  for(k = knz+1; k < blknnz; k++)
    if(cfound[k]){
      iwork[knz] = iwork[k];    /* block number (= new subscript) */
      ipos[knz] = ipos[k];      /* (old) position in yir */
      knz++;
    }
/* ------------------------------------------------------------
   Let blknnz be the number of nonzero trace entries
   ------------------------------------------------------------ */
  blknnz = knz;
  if(blknnz > 0){
/* ------------------------------------------------------------
   temporarily store nonzeros of tr-entries in fwork, then
   move other downzeros downwards to fill-up the gaps
   ------------------------------------------------------------ */
    for(knz = 0; knz < blknnz; knz++)
      fwork[knz] = ypr[ipos[knz]];
    for(knz = 1; knz < blknnz; knz++){
      inz = ipos[blknnz - knz - 1] + 1;     /* just beyond tr-entry */
      k = ipos[blknnz - knz] - inz;         /* amount before next tr-entry */
      if(k > 0){
        memmove(yir+inz+knz, yir+inz, k * sizeof(mwIndex));
        memmove(ypr+inz+knz, ypr+inz, k * sizeof(double));
      }
    }
    if(ipos[0] > 0){                   /* nonzeros before 1st nz-tr entry */
      memmove(yir+blknnz, yir, ipos[0] * sizeof(mwIndex));
      memmove(ypr+blknnz, ypr, ipos[0] * sizeof(double));
    }
/* ------------------------------------------------------------
   Re-insert the (saved) nonzero tr-entries at the start
   ------------------------------------------------------------ */
    memcpy(ypr, fwork, blknnz * sizeof(double));
    memcpy(yir, iwork, blknnz * sizeof(mwIndex));
    intadd(yir, blks[0], blknnz);
  }
}

/* ************************************************************
   PROCEDURE inverse of qreshape0, Transforms from [x1(lorN); x2-blocks]
     to pure block-wise [x[1], x[2],.. x[lorN]].
   INPUT
     blks - lorN array; blks[k]:blk[k+1]-1 are the subscripts of Lorentz
        block k=0:lorN-2.
     lorN - number of Lorentz blocks
   UPDATED
     y - entries in y(blk[0]:blks[lorN-1]) are affected (reshuffled)
   WORK
     fwork - length lorN-1 vector.
   ************************************************************ */
void qreshape1(double *y, const mwIndex *blks, const mwIndex lorN, double *fwork)
{
  mwIndex i,j,k;
  if(lorN <=1)
    return;             /* Nothing to do if only 1 block */
/* ------------------------------------------------------------
   Save y(2:lorN), which are the trace ("x1")-entries of block 2:lorN.
   ------------------------------------------------------------ */
  memcpy(fwork, y+blks[0] + 1, (lorN-1) * sizeof(double));
/* ------------------------------------------------------------
   For each block k = 0:lorN-2, shift y2[k] (lorN-1)-k positions upwards.
   ------------------------------------------------------------ */
  j = lorN;
  for(k = 1; k < lorN; k++){
    --j;                          /* lorN-1 >= j = lorN-k >= 1 */
    i = blks[k-1]+1;              /* new start of block y2[k-1] */
    memmove(y+i,y+i+j, (blks[k]-i) * sizeof(double));
  }
/* ------------------------------------------------------------
   Let y(iTr) = x(1:lorN), where iTr = cumsum([1 lorNL]) = blks.
   ------------------------------------------------------------ */
  for(k = 1; k < lorN; k++)
    y[blks[k]] = fwork[k-1];       /* y[blks[0]] is left unaffected */
}

/* ************************************************************
   PROCEDURE inverse of spqreshape0, Transforms from [x1(lorN); x2-blocks]
     to pure block-wise [x[1], x[2],.. x[lorN]].
   INPUT
     ynnz - nnz(y)
     blks - lorN+1 array; blks[0]:blks[1]-1 are the subscripts of the lorN
        trace values "y1[0:lorN-1]". Blks[k+1]:blk[k+2]-1 are the subscripts
        of Lorentz block y2[k], k = 0:lorN-2.
     lorN - number of Lorentz blocks
     iwsize - length of iwork, should be 2*(1+lorN).
   UPDATED
     yir, ypr - entries with subscripts in y(blk[0]:blks[lorN]) are
       affected (reshuffled)
   WORK
     cfound - length lorN character working array.
     iwork - length iwsize = 2 + lorN + MAX(floor(log(lorN+1)/log(2)), lorN-1),
        which is at most 2*(1+lorN).
     fwork - length lorN-1 vector.
   ************************************************************ */
/* still to be implemented */

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
 mwIndex i,j, ifirst, N,m, iwsize;
 coneK cK;
 mwIndex *iwork, *blks;
 double *fwork;
 bool *cwork;
 bool flag, bsparse;
 jcir y;

/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
 mxAssert(nrhs >= NPARIN, "qreshape requires more input arguments.");
 mxAssert(nlhs <= NPAROUT, "qreshape generates 1 output argument.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
 conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get inputs x, flag
   ------------------------------------------------------------ */
 N = mxGetM(X_IN);                    /* dimemsion */
 m = mxGetN(X_IN);                    /* number of cols */
 ifirst = 0;                          /* 1st Lorentz index */
 if(mxGetM(X_IN) != cK.qDim){
   mxAssert(mxGetM(X_IN) >= cK.lpN + cK.qDim, "x size mismatch");
   ifirst = cK.lpN;
 }
/* ------------------------------------------------------------
   Allocate output y(N,m)
   ------------------------------------------------------------ */
 Y_OUT = mxDuplicateArray(X_IN);                   /* Y = X_IN */
 if((bsparse = mxIsSparse(Y_OUT))){
   y.jc = mxGetJc(Y_OUT);
   y.ir = mxGetIr(Y_OUT);
 }
 y.pr = mxGetPr(Y_OUT);
 flag = (bool) mxGetScalar(FLAG_IN);
/* ------------------------------------------------------------
   Allocate working array iwork(2*(1+lorN)), fwork(lorN), blks(lorN),
   and cwork(lorN).
   ------------------------------------------------------------ */
 iwsize = 2 + 2 * cK.lorN;
 iwork = (mwIndex *) mxCalloc(iwsize, sizeof(mwIndex));
 fwork = (double *) mxCalloc(MAX(cK.lorN,1), sizeof(double));
 blks = (mwIndex *) mxCalloc(cK.lorN+1, sizeof(mwIndex));
 for(i = 1; i <= cK.lorN; i++)
   blks[i] = (mwIndex) cK.lorNL[i-1];           /* float to mwIndex */
 cwork = (bool *) mxCalloc(MAX(1,cK.lorN), sizeof(bool));
/* ------------------------------------------------------------
   Let iwork(0:lorN) = cumsum([ifirst, K.q(1:end)])
   ------------------------------------------------------------ */
 j = ifirst;
 blks[0] = j;
 for(i = 1; i <= cK.lorN; i++){
   j += blks[i];
   blks[i] = j;
 }
/* ------------------------------------------------------------
   The real job:
   ------------------------------------------------------------ */
 if(bsparse){
   mxAssert(flag != 1, "qreshape(x,1,K) cannot handle sparse x.");
#ifdef DO_ALL
   if(flag == 1)
     spqreshape1(y.pr, blks, cK.lorN, iwsize, iwork,fwork,cwork);
   else
#endif
     for(j = 0; j < m; j++)
       spqreshape0(y.ir+y.jc[j],y.pr+y.jc[j],y.jc[j+1]-y.jc[j],
                   blks, cK.lorN, iwsize, cwork, iwork,fwork);
 }
 else
   if(flag == 1)
     for(j = 0; j < m; j++)
       qreshape1(y.pr+j*N, blks, cK.lorN, fwork);
   else
     for(j = 0; j < m; j++)
       qreshape0(y.pr+j*N, blks, cK.lorN, fwork);
/* ------------------------------------------------------------
   RELEASE WORKING ARRAY
   ------------------------------------------------------------ */
 mxFree(iwork);
 mxFree(fwork);
 mxFree(cwork);
 mxFree(blks);
}
