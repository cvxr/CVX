/*
%                            Ad = Adendotd(dense, d, sparAd, Ablk, blkstart)
% ADENDOTD  Computes d[k]'*Aj[k] for Lorentz blocks that are to be factored
%  by dpr1fact.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function Ad = Adendotd(dense, d, sparAd, Ablk, blkstart)

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

#define AD_OUT plhs[0]
#define NPAROUT 1

#define DENSE_IN prhs[0]
#define D_IN prhs[1]
#define ADOTD_IN prhs[2]
#define ABLK_IN prhs[3]
#define BLKSTART_IN prhs[4]
#define NPARIN 5

/* ************************************************************
   PROCEDURE adendotd
   INPUT
     aden - dense.A(:,dense.l+1:end) sparse m x (nq+nden) matrix,
       with aden.jc[0] possibly nonzero.
     adotd - sparse m x nq array, has ai[k]'*d[k] for k in q, where
       the ai's are the sparse part of the A-matrix. We still need to
       add contribution from dense part, resulting in Ad.
     d1 - length |K.q| vector. We will use d1(q) entries.
     d2 - length firstQ+(sym(K.q)-|K.q|) vector. We use entries d2(dencols),
       where dencols >=firstQ.
     q - length nq array: dense lorentz blocks
     dencols - length nden array: dense lorentz norm-bound columns. These
       are global subscripts, at or beyond firstQ.
     blkend - length nq array, listing 1-beyond-last subscript of Lorentz
       norm bound blocks listed in q.
     nq - number of dense lorentz blocks
     nden - number of dense lorentz norm-bound columns
     fwork - length m vector.
   UPDATED
     ad - sparse m x nq. ad.ir and ad.jc are INPUTS, ad.pr is OUTPUT.
       On output, has (ai[k]+Adeni[k])'*d[k] for k in q.
   ************************************************************ */
void adendotd(jcir ad,jcir adotd,jcir aden,const double *d1,const double *d2,
              const mwIndex *q,const mwIndex *dencols,
              const mwIndex *blkend,const mwIndex nq,const mwIndex nden, double *fwork)
{
  mwIndex inz, i,j,k;
  const mwIndex *aden2jc;
  double dj;
/* ------------------------------------------------------------
   Initialize (Lorentz norm-bound part):
   1) aden2jc(0:nden) points to dense columns 
   2) j is next dense column to handle, inz point to next nonzero
   ------------------------------------------------------------ */
  j = 0;
  aden2jc = aden.jc + nq;   /* jump over Lorentz trace columns*/
  inz = aden2jc[j];
  for(k = 0; k < nq; k++){
/* ------------------------------------------------------------
   Set fwork = all-0;
   ------------------------------------------------------------ */
    for(i = ad.jc[k]; i < ad.jc[k+1]; i++)        /* fwork = all-0 */
      fwork[ad.ir[i]] = 0.0;
/* ------------------------------------------------------------
   Let fwork = adotd(:,k)    (Contribution from SPARSE part of A)
   ------------------------------------------------------------ */
    for(i = adotd.jc[k]; i < adotd.jc[k+1]; i++)
      fwork[adotd.ir[i]] = adotd.pr[i];
/* ------------------------------------------------------------
   Let fwork += d1(q(k)) * aden(:,k)   (Contribution Lorentz-trace)
   ------------------------------------------------------------ */
    dj = d1[q[k]];
    for(i = aden.jc[k]; i < aden.jc[k+1]; i++)
      fwork[aden.ir[i]] += dj * aden.pr[i];
/* ------------------------------------------------------------
   Add contribution of dense Lorentz-norm-bound columns, i.e.
   let fwork += sum_j{d2(dencols[j]) * Aden(:,j) | dencols[j]<blkend[k]}
   ------------------------------------------------------------ */
    for(; j < nden; j++){
      if((i = dencols[j]) >= blkend[k])
        break;                            /* Break if beyond block k */
      dj = d2[i];
      for(; inz < aden2jc[j+1]; inz++)
        fwork[aden.ir[inz]] += dj * aden.pr[inz];
    }
/* ------------------------------------------------------------
   Store ad(:,k) = fwork
   ------------------------------------------------------------ */
    for(i = ad.jc[k]; i < ad.jc[k+1]; i++)
      ad.pr[i] = fwork[ad.ir[i]];
  }
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  const mxArray *MY_FIELD;
  mwIndex i,j,firstQ, m,nden, nl, nq, lorN;
  mwIndex *q, *dencols, *blkend;
  const double *d1, *d2, *qPr, *dencolsPr, *blkstartPr;
  double *fwork;
  jcir ad, aden,adotd;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "adendotd requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "adendotd produces less output arguments");
/* ------------------------------------------------------------
   DISASSEMBLE dense structure: dense.{cols,l,q,A}
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(DENSE_IN),"dense should be a structure.");
  MY_FIELD = mxGetField(DENSE_IN,(mwIndex)0,"l");               /* dense.l */
  mxAssert( MY_FIELD != NULL, "Missing field dense.l.");
  nl = (mwIndex) mxGetScalar(MY_FIELD);                             /* double to mwIndex */
  MY_FIELD = mxGetField(DENSE_IN,(mwIndex)0,"q");                    /* dense.q */
  mxAssert( MY_FIELD != NULL, "Missing field dense.q.");
  nq = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
  qPr = mxGetPr(MY_FIELD);
  MY_FIELD = mxGetField(DENSE_IN,(mwIndex)0,"cols");                 /* dense.cols */
  mxAssert( MY_FIELD != NULL, "Missing field dense.cols.");
  nden = mxGetM(MY_FIELD) * mxGetN(MY_FIELD) - nl - nq;
  mxAssert(nden >= 0, "dense.q size mismatch.");
  dencolsPr = mxGetPr(MY_FIELD) + nl + nq;               /* Skip LP and Q-tr*/
  MY_FIELD = mxGetField(DENSE_IN,(mwIndex)0,"A");                     /* dense.A */
  mxAssert( MY_FIELD != NULL, "Missing field dense.A.");
  mxAssert(mxIsSparse(MY_FIELD), "dense.A must be sparse");
  m = mxGetM(MY_FIELD); 
  mxAssert(mxGetN(MY_FIELD) - nl == nq + nden, "dense.A size mismatch");
  aden.jc = mxGetJc(MY_FIELD) + nl;                         /* Skip LP part */
  aden.ir = mxGetIr(MY_FIELD);
  aden.pr = mxGetPr(MY_FIELD);
/* ------------------------------------------------------------
   DISASSEMBLE d structure: d.{q1,q2}
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(D_IN), "d should be a structure.");
  MY_FIELD = mxGetField(D_IN,(mwIndex)0,"q1");         /* d.q1 */
  mxAssert( MY_FIELD != NULL, "Missing field d.q1.");
  lorN = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
  d1 = mxGetPr(MY_FIELD);
  MY_FIELD = mxGetField(D_IN,(mwIndex)0,"q2");        /* d.q2 */
  mxAssert( MY_FIELD != NULL, "Missing field d.q2.");
  d2 = mxGetPr(MY_FIELD);
/* ------------------------------------------------------------
   Get inputs adotd (contains Ad from sparse A in dense.qs blocks),
   blkstart (partitions d2 into Lorentz norm-bound blocks)
   ------------------------------------------------------------ */
  mxAssert(mxIsSparse(ADOTD_IN), "sparAD must be sparse");      /* adotd */
  mxAssert((m == mxGetM(ADOTD_IN) || nq <= 0) && nq == mxGetN(ADOTD_IN), "Size mismatch sparAD");
  adotd.jc = mxGetJc(ADOTD_IN);
  adotd.ir = mxGetIr(ADOTD_IN);
  adotd.pr = mxGetPr(ADOTD_IN);
  blkstartPr = mxGetPr(BLKSTART_IN);                  /* blkstart */
  mxAssert(lorN +1 == mxGetM(BLKSTART_IN) * mxGetN(BLKSTART_IN), "blkstart size mismatch");
/* ------------------------------------------------------------
   Create working arrays q(nq), dencols(nden), fwork(m),
   blkend(nq)
   ------------------------------------------------------------ */
  q    = (mwIndex *) mxCalloc(MAX(1,nq), sizeof(mwIndex));
  dencols = (mwIndex *) mxCalloc(MAX(1,nden), sizeof(mwIndex));
  blkend = (mwIndex *) mxCalloc(MAX(1,nq), sizeof(mwIndex));
  fwork = (double *) mxCalloc(MAX(m,1), sizeof(double));
/* ------------------------------------------------------------
   Convert to integer C-style; dencols, q, blkstart(q+1)
   ------------------------------------------------------------ */
  for(i = 0; i < nden; i++){
    j = (mwIndex) dencolsPr[i];
    mxAssert(j>0,"");
    dencols[i] = --j;
  }
  for(i = 0; i < nq; i++){
    j = (mwIndex) qPr[i];
    mxAssert(j>0,"");
    q[i] = --j;
  }
/* ------------------------------------------------------------
   Let firstQ point to subscript of 1st Lorentz norm-bound variable
   ------------------------------------------------------------ */
  firstQ = (mwIndex) blkstartPr[0];            /* double to mwIndex */
  mxAssert(firstQ>0,"");  
  --firstQ;                          /* Fortran to C */
  for(i = 0; i < nq; i++){
    j = (mwIndex) blkstartPr[q[i] + 1];        /* F-double to C-mwIndex */
    mxAssert(j>0,"");
    blkend[i] = --j;
  }
/* ------------------------------------------------------------
   Create output: Ad = Ablk
   ------------------------------------------------------------ */
  AD_OUT = mxDuplicateArray(ABLK_IN);              /* Ad = Ablk */
  ad.jc = mxGetJc(AD_OUT);
  ad.ir = mxGetIr(AD_OUT);
  ad.pr = mxGetPr(AD_OUT);
/* ------------------------------------------------------------
   The real job is done here:
   ------------------------------------------------------------ */
  adendotd(ad,adotd,aden,d1,d2 - firstQ,q,dencols,blkend,nq,nden, fwork);
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(fwork);
  mxFree(dencols);
  mxFree(q);
}
