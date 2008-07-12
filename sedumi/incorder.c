/* ************************************************************
%                            [perm, dz] = incorder(At [,Ajc1,ifirst])
% INCORDER
% perm sorts the columns of At greedily, by iteratively picking
%   the 1st unprocessed column with the least number of nonzero
%   subscripts THAT ARE NOT YET COVERED (hence incremental) by
%   the previously processed columns.
% dz has the corresponding incremental sparsity structure, i.e.
%   each column lists only the ADDITIONAL subscripts w.r.t. the
%   previous ones.
%
% SEE ALSO getada3, dpr1fact.
% ******************** INTERNAL FUNCTION OF SEDUMI ********************

function [perm, dz] = incorder(At,Ajc1,ifirst)
   
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

#define PERM_OUT myplhs[0]    /* incremental ordering */
#define DZ_OUT myplhs[1]      /* incremental nonzero structure */
#define NPAROUT 2

#define AT_IN prhs[0]
#define NPARINMIN 1
#define AJC1_IN prhs[1]
#define IFIRST_IN prhs[2]
#define NPARIN 3

/* ************************************************************
   PROCEDURE spPartTransp - Let (Ajc, Air) = At(first:lenfull-1,:)'.
   INPUT
     Atir - length Atjc2[m-1] subscript array.
     Atjc1 - length m array, pointing to 1st index in At(:,j) that is
       in range first:lenfull-1.
     Atjc2 - length m array, pointing to column end.
     first, lenfull - index range we're interested in. (lenud:=lenfull-first)
     m - number of constraints (rows in At).
   OUTPUT
     Ajc - length 1+lenud mwIndex-array, column pointers
       of A := At(first:lenfull-1,:)'.
     Air - subscript array of output matrix A. length is:
       sum(Atjc(2:m) - Atjcs(1:m)).
   WORK
     iwork - length lenud (=lenfull-first) integer array.
   ************************************************************ */
void spPartTransp(mwIndex *Air, mwIndex *Ajc,
                  const mwIndex *Atir, const mwIndex *Atjc1, const mwIndex *Atjc2,
                  const mwIndex first, const mwIndex lenfull, const mwIndex m,
                  mwIndex *iwork)
{
  mwIndex i,j,inz;
  mwIndex *inxnz;
/* ------------------------------------------------------------
   INIT: Make indices Ajc[first:lenfull] valid. Then set Ajc(:)=all-0.
   ------------------------------------------------------------ */
  Ajc -= first;
  for(i = first; i <= lenfull; i++)
    Ajc[i] = 0;
/* ------------------------------------------------------------
   For each column j=0:m-1, in each nz PSD entry Atj(i): ++inxnz[i],
   where inxnz := Ajc+1. (we use Ajc[first+1:lenfull])
   ------------------------------------------------------------ */
  inxnz = Ajc + 1;
  for(j = 0; j < m; j++)
    for(inz = Atjc1[j]; inz < Atjc2[j]; inz++)
      ++inxnz[Atir[inz]];
/* ------------------------------------------------------------
   cumsum(inxnz). Note that Ajc[0] =0 already.
   ------------------------------------------------------------ */
  for(inz = 0, i = first; i < lenfull; i++)
    inxnz[i] += inxnz[i-1];
/* ------------------------------------------------------------
   Write the subscripts of A:= At(first:lenfull-1,:)'.
   ------------------------------------------------------------ */
  memcpy(iwork, Ajc + first, (lenfull - first) * sizeof(mwIndex));
  iwork -= first;          /* points to next avl. entry in A(:,i) */
  for(j = 0; j < m; j++)
    for(inz = Atjc1[j]; inz < Atjc2[j]; inz++){
      i = Atir[inz];
      Air[iwork[i]++] = j;         /* as (j,i) entry */
    }
}

/* ************************************************************
   PROCEDURE incorder - Greedy ordering of PSD-part of columns in At,
       starting with least number of (PSD)-nonzeros. Good fot getada3.
       Produces N x m sparse matrix dz, having at most lenud nonzeros.
       The last columns correspond to the "deg" columns, and are not
       reordered. The first m-ndeg are reordered, see "perm".
   INPUT
     Atjc1 - length m, start of At(first:lenful,j) for each j=1:m.
     Atjc2 - column end pointers, length m.
     Atir - row-subscripts of At.
     Ajc, Air - sparsity structure of A := At(first:lenful,:)', is m x lenud.
     m    - number of columns in At.
     first - first PSD subscipt.
   OUTPUT
     perm - length m array. perm(1:m-ndeg) is the greedy ordering,
       starting from sparsest. perm(m-ndeg:m) = deg.
     dzjc - length m+1, column pointers into dzir.
     dzir - length dzjc[m] <= lenud. jth column lists additional
       subscripts to dz(:,0:j-1).
   WORK
     iwork - (length m) Remaining length of PSD part
       in each column j=1:m.
     discard - length lenud of Booleans. 1 means index already in dz.
   ************************************************************ */
void incorder(mwIndex *perm, mwIndex *dzjc, mwIndex *dzir,
              const mwIndex *Atjc1,const mwIndex *Atjc2,const mwIndex *Atir,
              const mwIndex *Ajc,const mwIndex *Air,
              const mwIndex m, const mwIndex first, const mwIndex lenud,
              mwIndex *iwork, char *discard)
{
  mwIndex kmin,lenmin,i,j,k, inz,jnz, permk;
/* ------------------------------------------------------------
   initialize: dzjc[0]=0, nperm = m - ndeg. Let Ajc-=first, so
   that Ajc[first:lenfull] is valid. Similar for discard.
   Initialize discard to all-0 ("False").
   ------------------------------------------------------------ */
  dzjc[0] = 0;
  Ajc -= first;
  memset(discard, '\0', lenud);            /* all-0 */
  discard -= first;
/* ------------------------------------------------------------
   Let perm(1:m) be 1:m
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++)
    perm[i] = i;
/* ------------------------------------------------------------
   Set iwork = Atjc2(1:m)-Atjc1(1:m).
   The number of (not yet discarded) nz-PSD indices per constraint.
   ------------------------------------------------------------ */
  for(k = 0; k < m; k++){
    iwork[k] = Atjc2[k] - Atjc1[k];
  }
/* ------------------------------------------------------------
   In iterate k=0:m-1, pivot on constraint with length iwork[k] minimal.
   ------------------------------------------------------------ */
  for(k = 0; k < m; k++){
    kmin = k;
    lenmin = iwork[perm[k]];
    for(j = k+1; j < m; j++)
      if(iwork[perm[j]] < lenmin){
        lenmin = iwork[perm[j]];
        kmin = j;
      }
    mxAssert(lenmin >= 0,"");
    permk = perm[kmin];                  /* make pivot in perm */
    perm[kmin] = perm[k];
    perm[k] = permk;
/* ------------------------------------------------------------
   Write the (additional) subscripts of the pivot column into
   k-th column of dz, i.e. dz(:,k) = At(:,permk).
   ------------------------------------------------------------ */
    jnz = dzjc[k];
    for(inz = Atjc1[permk]; inz < Atjc2[permk]; inz++){
      i = Atir[inz];
      if(!discard[i]){
        discard[i] = 1;
        dzir[jnz++] = i;
      }
    }
    mxAssert(jnz == dzjc[k] + lenmin,"");
/* ------------------------------------------------------------
   Discard dz(:,k)'s subscripts: adjust the constraint lengths
   where applicable.
   ------------------------------------------------------------ */
    dzjc[k+1] = jnz;
    for(jnz = dzjc[k]; jnz < dzjc[k+1]; jnz++){
      i = dzir[jnz];                         /* i is discarded subscript */
      for(inz = Ajc[i]; inz < Ajc[i+1]; inz++){
        j = Air[inz];
        iwork[j]--;                   /* subscript discard here */
      }
    }
  }
}


/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
      [perm,dz] = incorder(At,Ajc1,ifirst)   
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  mxArray *myplhs[NPAROUT];
  mwIndex i, m, firstPSD, lenud, iwsiz, lenfull, maxnnz;
  mwIndex *iwork, *Atjc1, *Ajc, *Air, *perm;
  const mwIndex *Atjc2, *Atir;
  double *permPr;
  const double *Ajc1Pr;
  char *cwork;
  jcir dz;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARINMIN, "incorder requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "incorder produces less output arguments.");
/* --------------------------------------------------
   GET STATISTICS:
   -------------------------------------------------- */
  if(nrhs == NPARIN){
    firstPSD = (mwIndex) mxGetScalar(IFIRST_IN);        /* double to mwIndex */
    mxAssert(firstPSD>0,"");
    --firstPSD;                               /* Fortran to C */
  }
  else
    firstPSD = 0;                           /* default: use all subscripts */
  lenfull = mxGetM(AT_IN);
  lenud = lenfull - firstPSD;
/* --------------------------------------------------
   Check size At, and get At.
   -------------------------------------------------- */
  mxAssert(mxIsSparse(AT_IN), "At must be a sparse matrix.");
  m = mxGetN(AT_IN);
  Atjc2 = mxGetJc(AT_IN)+1;        /* points to end of constraint */
  Atir = mxGetIr(AT_IN);           /* subscripts */
/* ------------------------------------------------------------
   ALLOCATE WORKING arrays:
   mwIndex Atjc1(m),  Ajc(1+lenud), Air(maxnnz), perm(m),
     iwork(iwsiz). iwsiz = MAX(lenud,1+m)
   char cwork(lenud).
   ------------------------------------------------------------ */
  Atjc1 = (mwIndex *) mxCalloc(MAX(m,1), sizeof(mwIndex));
  Ajc = (mwIndex *) mxCalloc(1+lenud, sizeof(mwIndex));
  perm = (mwIndex *) mxCalloc(MAX(1,m), sizeof(mwIndex));
  iwsiz = MAX(lenud, 1 + m);
  iwork = (mwIndex *) mxCalloc(iwsiz, sizeof(mwIndex)); /* iwork(iwsiz) */
  cwork = (char *) mxCalloc(MAX(1,lenud), sizeof(char));
/* ------------------------------------------------------------
   Get input AJc1
   ------------------------------------------------------------ */
  if(nrhs == NPARIN){
    Ajc1Pr = mxGetPr(AJC1_IN);
    mxAssert(mxGetM(AJC1_IN) * mxGetN(AJC1_IN) >= m, "Ajc1 size mismatch");
/* ------------------------------------------------------------
   Double to mwIndex: Atjc1
   ------------------------------------------------------------ */
    for(i = 0; i < m; i++)
      Atjc1[i] = (mwIndex) Ajc1Pr[i];           /* double to mwIndex */
/* ------------------------------------------------------------
   Let maxnnz = number of PSD nonzeros = sum(Atjc2-Atjc1)
   ------------------------------------------------------------ */
  maxnnz = 0;
  for(i = 0; i < m; i++)
    maxnnz += Atjc2[i] - Atjc1[i];
  }
  else{
    memcpy(Atjc1, mxGetJc(AT_IN), m * sizeof(mwIndex)); /* default column start */
    maxnnz = Atjc2[m-1];                    /* maxnnz = nnz(At) */
  }
/* ------------------------------------------------------------
   ALLOCATE WORKING array mwIndex Air(maxnnz)
   ------------------------------------------------------------ */
  Air = (mwIndex *) mxCalloc(MAX(1,maxnnz), sizeof(mwIndex));
/* ------------------------------------------------------------
   CREATE OUTPUT ARRAYS PERM(m) and DZ=sparse(lenfull,m,lenud).
   ------------------------------------------------------------ */
  PERM_OUT = mxCreateDoubleMatrix(m,(mwSize)1,mxREAL);
  permPr = mxGetPr(PERM_OUT);
  DZ_OUT = mxCreateSparse(lenfull,m,lenud,mxREAL);
  dz.jc = mxGetJc(DZ_OUT);
  dz.ir = mxGetIr(DZ_OUT);
/* ------------------------------------------------------------
   Let (Ajc,Air) := At(first:end,:)', the transpose of PSD-part.
   Uses iwork(lenud)
   ------------------------------------------------------------ */
  spPartTransp(Air,Ajc, Atir,Atjc1,Atjc2, firstPSD,lenfull, m,iwork);
/* ------------------------------------------------------------
   The main job: greedy order of columns of At
   ------------------------------------------------------------ */
  incorder(perm, dz.jc,dz.ir, Atjc1,Atjc2,Atir, Ajc,Air,
           m, firstPSD,lenud, iwork, cwork);    /* uses iwork(m+1) */
/* ------------------------------------------------------------
   REALLOC (shrink) dz to dz.jc[m] nonzeros.
   ------------------------------------------------------------ */
  mxAssert(dz.jc[m] <= lenud,"");
  maxnnz = MAX(1,dz.jc[m]);                     /* avoid realloc to 0 */
  if((dz.ir = (mwIndex *) mxRealloc(dz.ir, maxnnz * sizeof(mwIndex))) == NULL)
    mexErrMsgTxt("Memory allocation error");
  if((dz.pr = (double *) mxRealloc(mxGetPr(DZ_OUT), maxnnz*sizeof(double)))
     == NULL)
    mexErrMsgTxt("Memory allocation error");
  mxSetPr(DZ_OUT,dz.pr);
  mxSetIr(DZ_OUT,dz.ir);
  mxSetNzmax(DZ_OUT,maxnnz);
  for(i = 0; i < maxnnz; i++)
    dz.pr[i] = 1.0;
/* ------------------------------------------------------------
   Convert C-mwIndex to Fortran-double
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++)
    permPr[i] = perm[i] + 1.0;
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(cwork);
  mxFree(iwork);
  mxFree(perm);
  mxFree(Air);
  mxFree(Ajc);
  mxFree(Atjc1);
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  i = MAX(nlhs, 1);
  memcpy(plhs,myplhs, i * sizeof(mxArray *));
  for(; i < NPAROUT; i++)
    mxDestroyArray(myplhs[i]);
}
