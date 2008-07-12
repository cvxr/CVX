/* ************************************************************
%                                       ADA = getada2(ADA, DAt,Aord, K)
% GETADA2  Compute ADA += DAt.q'*DAt.q
%   IMPORTANT: Updated ADA only on triu(ADA(Aord.qperm,Aord.qperm)).
%     Remaining entries are not affected.
%
% SEE ALSO sedumi, getada1, getada3
% ******************** INTERNAL FUNCTION OF SEDUMI ********************

function ADA = getada2(ADA, DAt,Aord, K)

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
#include <math.h>
#include "mex.h"
#include "blksdp.h"

#define ADA_OUT plhs[0]
#define NPAROUT 1

#define ADA_IN prhs[0]       /* sparsity struct ADA */
#define DAT_IN prhs[1]       /* structure with DAt.q (=ddota) */
#define AORD_IN prhs[2]
#define K_IN  prhs[3]
#define NPARIN 4

/* ************************************************************
   PROCEDURE: getada2 - Let ADA += ddota'*ddota.
   INPUT
     ada.{jc,ir} - sparsity structure of ada.
     ddota - sparse lorN x m matrix.
     perm, invperm - length(m) array, ordering in which ADA should be computed,
       and its inverse. We compute in order triu(ADA(perm,perm)), but store
       at original places. OPTIMAL PERM: sort(sum(spones(ddota))), i.e. start
       with sparsest.
     m  - order of ADA, number of constraints.
     lorN - length(K.q), number of Lorentz blocks.
   UPDATED
     ada.pr - ada(i,j) += ddotai'*ddotaj. ONLY triu(ADA(perm,perm)) is
        updated. (So caller typically should symmetrize afterwards.)
   WORKING ARRAYS
     ddotaj - work vector, size lorN.
   ************************************************************ */
void getada2(jcir ada, jcir ddota, const mwIndex *perm, const mwIndex *invperm,
             const mwIndex m, const mwIndex lorN,   double *ddotaj)
{
  mwIndex i,j, knz,inz, permj;
  double adaij;

/* ------------------------------------------------------------
   Init ddotaj = all-0 (for Lorentz)
   ------------------------------------------------------------ */
  fzeros(ddotaj, lorN);
/* ============================================================
   MAIN getada LOOP: loop over nodes perm(0:m-1)
   ============================================================ */
  for(j = 0; j < m; j++){
    permj = perm[j];
    if(ddota.jc[permj] < ddota.jc[permj+1]){      /* Only work if nonempty */
/* ------------------------------------------------------------
   Let ddotaj = ddota(:,j) in full
   ------------------------------------------------------------ */
      for(i = ddota.jc[permj]; i < ddota.jc[permj+1]; i++)
        ddotaj[ddota.ir[i]] = ddota.pr[i];
/* ------------------------------------------------------------
   For all i with invpermi < j:
   ada_ij += ddota_i'*ddotaj.
   ------------------------------------------------------------ */
      for(inz = ada.jc[permj]; inz < ada.jc[permj+1]; inz++){
        i = ada.ir[inz];
        if(invperm[i] <= j){
          adaij = ada.pr[inz];
          if(invperm[i] < j)
            for(knz = ddota.jc[i]; knz < ddota.jc[i+1]; knz++)
              adaij +=  ddota.pr[knz] * ddotaj[ddota.ir[knz]];
          else                         /* diag entry: += ||ddota(:,permj)||^2 */
            adaij += realssqr(ddota.pr + ddota.jc[i], ddota.jc[i+1]-ddota.jc[i]);
          ada.pr[inz] = adaij;
        }
      }
/* ------------------------------------------------------------
   Re-initialize ddotaj = 0.
   ------------------------------------------------------------ */
      for(i = ddota.jc[permj]; i < ddota.jc[permj+1]; i++)      /* Lorentz */
        ddotaj[ddota.ir[i]] = 0.0;
    }
  } /* j = 0:m-1 */
}


/* ============================================================
   MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  coneK cK;
  const mxArray *MY_FIELD;
  mwIndex m, i, j;
  const double *permPr;
  double *fwork;
  mwIndex *iwork, *perm, *invperm;
  jcir ada, ddota;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "getADA requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "getADA produces less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
  m = mxGetM(ADA_IN);
/* ------------------------------------------------------------
   Allocate output matrix ADA with sparsity structure of ADA_IN,
   and initialize as a copy of ADA_IN.
   ------------------------------------------------------------ */
  mxAssert(mxGetN(ADA_IN) == m, "Size mismatch ADA.");
  mxAssert(mxIsSparse(ADA_IN), "ADA should be sparse.");
  ADA_OUT = mxDuplicateArray(ADA_IN);              /* ADA = ADA_IN */
  if(cK.lorN <= 0)                          /* READY if no LORENTZ blocks !*/
    return;
  ada.jc = mxGetJc(ADA_OUT);
  ada.ir = mxGetIr(ADA_OUT);
  ada.pr = mxGetPr(ADA_OUT);
/* ------------------------------------------------------------
   DISASSEMBLE DAt structure: DAt.q
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(DAT_IN), "DAt should be a structure.");
  MY_FIELD = mxGetField(DAT_IN,(mwIndex)0,"q");       /* DAt.q */
  mxAssert( MY_FIELD != NULL, "Missing field DAt.q.");
  mxAssert(mxGetM(MY_FIELD) == cK.lorN && mxGetN(MY_FIELD) == m, "Size mismatch DAt.q");
  mxAssert(mxIsSparse(MY_FIELD), "DAt.q should be sparse.");
  ddota.jc = mxGetJc(MY_FIELD);
  ddota.ir = mxGetIr(MY_FIELD);
  ddota.pr = mxGetPr(MY_FIELD);
/* ------------------------------------------------------------
   DISASSEMBLE Aord structure: Aord.qperm
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(AORD_IN), "Aord should be a structure.");
  MY_FIELD = mxGetField(AORD_IN,(mwIndex)0,"qperm");        /* Aord.qperm */
  mxAssert( MY_FIELD != NULL, "Missing field Aord.qperm.");
  mxAssert(mxGetM(MY_FIELD) * mxGetN(MY_FIELD) == m, "Size mismatch Aord.qperm.");
  permPr = mxGetPr(MY_FIELD);
/* ------------------------------------------------------------
   Only work to do if ~isempty(ddota):
   ------------------------------------------------------------ */
  if(ddota.jc[m] > 0){
/* ------------------------------------------------------------
   ALLOCATE working arrays:
   iwork(2*m) = [perm(m), invperm(m)].
   fwork[lorN]
   ------------------------------------------------------------ */
    iwork = (mwIndex *) mxCalloc(MAX(2 * m,1), sizeof(mwIndex));
    perm = iwork;
    invperm = perm + m;
    fwork  = (double *) mxCalloc(MAX(cK.lorN,1), sizeof(double));
/* ------------------------------------------------------------
   perm to integer C-style
   ------------------------------------------------------------ */
    for(i = 0; i < m; i++){
      j = (mwIndex) permPr[i];
      mxAssert(j>0,"");
      perm[i] = --j;
    }
/* ------------------------------------------------------------
   Let invperm(perm) = 0:m-1.
   ------------------------------------------------------------ */
    for(i = 0; i < m; i++)
      invperm[perm[i]] = i;
/* ------------------------------------------------------------
   ACTUAL COMPUTATION: ADA += DAt.q'*DAt.q.
   ------------------------------------------------------------ */
    getada2(ada, ddota, perm, invperm, m, cK.lorN, fwork);
/* ------------------------------------------------------------
   RELEASE WORKING ARRAYS.
   ------------------------------------------------------------ */
    mxFree(fwork);
    mxFree(iwork);
  } /* !isempty(ddota) */
}
