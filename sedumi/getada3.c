/* ************************************************************
%                            [ADA,absd] = getada3(ADA, A,Ajc1,Aord, udsqr,K)
% GETADA3  Compute ADA(i,j) = (D(d^2)*A.t(:,i))' *A.t(:,j),
%   and exploit sparsity as much as possible.
%   absd - length m output vector, containing
%     absd(i) = abs((D(d^2)*A.t(:,i))' *abs(A.t(:,i)).
%     Hence, diag(ADA)./absd gives a measure of cancelation (in [0,1]).
%
% SEE ALSO sedumi, getada1, getada2
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
  
   function [ADA,absd] = getada3(ADA, A,Ajc1,Aord, udsqr,K)

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

#define ADA_OUT myplhs[0]
#define ABSD_OUT myplhs[1]
#define NPAROUT 2

#define ADA_IN prhs[0]       /* sparsity struct ADA */
#define AT_IN prhs[1]        /* N x m sparse At */
#define AJC1_IN prhs[2]      /* start of PSD blocks in At */
#define AORD_IN prhs[3]
#define UDSQR_IN prhs[4]     /* Scale data */
#define K_IN  prhs[5]
#define NPARIN 6

/* ========================= G E N E R A L ========================= */

void spzeros(double *x,const mwIndex *xir,const mwIndex n)
{
  mwIndex i;
  for(i = 0; i < n; i++)
    x[xir[i]] = 0.0;
}

/* ************************************************************
   PROCEDURE exmerge - mergeing 2 exclusive, increasing integer arrays.
   INPUT
     x - length nx array, increasing entries.
     y - length ny array, its entries are increasing, and do not occur in x.
     nx,ny - order of x and y.
   OUTPUT
     z - length nx+ny vector
   WORK
     cwork - length ny
     iwork - length ny+2+floor(log_2(1+ny)).
   ************************************************************ */
void exmerge(mwIndex *x, const mwIndex *y, const mwIndex nx, const mwIndex ny,
             const mwIndex iwsize, bool *cwork, mwIndex *ipos)
{
  mwIndex i,j,inz;
/* ------------------------------------------------------------
   Search all "insertion" positions of y-entries in list x
   ------------------------------------------------------------ */
  intmbsearch(ipos,cwork, x,nx,y,ny, ipos+ny+2,iwsize-(ny+2));
/* ------------------------------------------------------------
   Shift x-entries down to make room for y-entries
   ------------------------------------------------------------ */
  for(i = ny; i > 0; i--){         /* shift i entries down */
    inz = ipos[i];
    if((j = ipos[i+1]-inz) > 0)
      memmove(x+inz+i,x+inz,j*sizeof(mwIndex));
  }
/* ------------------------------------------------------------
   Insert y-entries
   ------------------------------------------------------------ */
  ++ipos;
  for(i = 0; i < ny; i++)
    x[ipos[i]+i] = y[i];          /* add i, number prev. inserted y's */
#ifndef NDEBUG
  for(i = 1; i < nx+ny; i++)
    mxAssert(x[i] > x[i-1],"");
#endif
}

/* ************************************************************
   PROCEDURE cpspdiag -- Let d := diag(X), where X is sparse m x m.
   INPUT
     m - number of columns in ada
     x - contains COMPLETE SYMMETRIC sparsity structure,
   OUTPUT
     d - length m vector, diag(X).
   ************************************************************ */
void cpspdiag(double *d, const jcir x, const mwIndex m)
{
  mwIndex j, inz;
  const mwIndex *diagFound;
  
/* ------------------------------------------------------------
   For each column j: let dj = x(j,j)
   ------------------------------------------------------------ */
  for(j = 0; j < m; j++){
    if((diagFound = ibsearch(&j,x.ir+x.jc[j],x.jc[j+1]-x.jc[j])) != NULL){
      inz = diagFound - x.ir;                 /* convert pointer to index */
      d[j] = x.pr[inz];                       /* diag entry found */
    }
    else
      d[j] = 0.0;                   /* no diag entry -> d[j] = 0.0 */
  }
}

/* ************************************************************
   PROCEDURE spmakesym -- Let X := X+X'-diag(diag(X)), with X a
      real sparse square matrix.
   INPUT
     m - number of columns in ada
   UPDATED
     x - on input, contains COMPLETE SYMMETRIC sparsity structure,
         and the values (possibly 0's on some locations). On return,
         X = X+X'-diag(diag(X)), i.e. X is symmetrized, and the off-
         diagonal entries are "DUPLICATED" (allowing part in xij and xji).
   WORK
     iwork - length m array of integers. Points to "below row j"
       part of columns (trilstart). (initial contents irrelevant)
   ************************************************************ */
void spmakesym(jcir x, const mwIndex m, mwIndex *iwork)
{
  mwIndex i, j, inz, jend;
  double xij;
  
/* ------------------------------------------------------------
   Initialize: let iwork(0:m-1) = ada.jc(0:m-1)
   ------------------------------------------------------------ */
  memcpy(iwork, x.jc, m * sizeof(mwIndex));   /* don't copy x.jc[m] */
/* ------------------------------------------------------------
   For each column j:   for each index i > j:
   let xij = x(i,j) + x(j,i).
   Let iwork point to next nonzero in col i
   Guard: x.ir[iwork(i)] >= j for all i >= j.
   ------------------------------------------------------------ */
  for(j = 0; j < m; j++){
    jend = x.jc[j+1];                        /* Let [inz,jend) be tril(:,j) */
    inz = iwork[j];
    if(inz < jend){
      if(x.ir[inz] == j)                         /* skip diagonal entry */
        inz++;
      while(inz < jend){                         /* off-diagonal entries */
        i = x.ir[inz];
        xij = x.pr[iwork[i]] + x.pr[inz];  /* xji + xij */
        x.pr[iwork[i]++] = xij;
        x.pr[inz++] = xij;
      }
    }
  }
}

/* ************************************************************
   PROCEDURE dzblkpartit
   INPUT
     dzir - length dznnz array, NOT sorted
     xblk - length max(dzir) array, xblk(dzir) maps into 0:nblk-1.
     dznnz - order of dzir
     nblk - 1+max(xblk), number of blocks
   OUTPUT
     dzjc - length nblk+1 array. Has blockstarts so that all subscripts
       in dzir fit in the resulting partition.
   ************************************************************ */
void dzblkpartit(mwIndex *dzjc, const mwIndex *dzir, const mwIndex *xblk,
                 const mwIndex dznnz, const mwIndex nblk)
{
  mwIndex i,j;
/* ------------------------------------------------------------
   Init dzjc = all-0
   ------------------------------------------------------------ */
  for(i = 0; i <= nblk; i++)
    dzjc[i] = 0;
/* ------------------------------------------------------------
   Pre-partition dzir(dznnz) space into blocks
   ------------------------------------------------------------ */
  ++dzjc;
  for(i = 0; i < dznnz; i++)
    dzjc[xblk[dzir[i]]]++;      /* accumulate nnz */
  j = 0;
  for(i = 0; i < nblk; i++){    /* cumsum */
    j += dzjc[i];
    dzjc[i] = j;
  }
  mxAssert(dzjc[nblk-1] == dznnz,"");
}    


/* ************************************************************
   PROCEDURE: getada3
   INPUT
     ada.{jc,ir} - sparsity structure of ada.
     At - sparse N x m matrix.
     udsqr - lenud vector containing D, D(ud.perm,ud.perm) = Ud'*Ud.
     Ajc1 - m mwIndex array, Ajc1(:,1) points to start of PSD nz's in At.
     dzjc - psdN+1, partition of dz rowsubscipts into PSD blocks.
     dzstructjc, dzstructir - sparse N x m matrix, giving NEW PSD-nonzero
       positions of At(:,perm(j)).
     blkstart - length(K.s): starting indices of PSD  blocks
     xblk - length psdDim array, with k = xblk(i-blkstart[0]) iff
       blkstart[k] <= i < blkstart[k+1], k=0:nblk-1.
       psdDim:=blkstart[end]-blkstart[0].
     psdNL - K.s in integer
     perm, invperm - length(m) array, ordering in which ADA should be computed,
       and its inverse. We compute in order triu(ADA(perm,perm)), but store
       at original places. OPTIMAL order: start with sparsest dzstruct.
     m  - order of ADA, number of constraints.
     lenud - blkstart[end] - blkstart[0], PSD dimension.
     pcK - pointer to cone K structure.
     rpsdN - sum(K.s(i) | i is real sym PSD block).
   UPDATED
     ada.pr - ada(i,j) += ai'*D(d^2; PSD)*aj. PSD-part. Only entries
       in triu(ADA(perm,perm) are affected.
   OUTPUT
     absd - length(m) vector, contains abs(aj)'*abs(D(d^2)*aj).
   WORKING ARRAYS
     fwork - work vector, size
         fwsiz = lenud + 2 * max(rMaxn^2, 2*hMaxn^2).
     iwork - integer work array, size 
         iwsiz = 2*psdN + dznnz + max(srchsize, max(nk(PSD))).
	 where dznnz = dzstructjc[m] and
         srchsize = 2+maxadd+floor(log(1+maxadd))
     cwork - maxadd, where maxadd = max(dzstructjc(i+1)-dzstructjc(i))
   ************************************************************ */
void getada3(jcir ada, double *absd, jcir At, const double *udsqr,
             const mwIndex *Ajc1, const mwIndex *dzjc,
             const mwIndex *dzstructjc, const mwIndex *dzstructir,
             const mwIndex *blkstart, const mwIndex *xblk, const mwIndex *psdNL,
             const mwIndex *perm, const mwIndex *invperm,
             const mwIndex m, const mwIndex lenud, const coneK *pcK,
             double *fwork, mwIndex fwsiz, mwIndex *iwork, mwIndex iwsiz,
             bool *cwork)
{
  mwIndex i,j,k, knz,inz, dznnz, permj, rsdpN, nblk, nnzbj;
  double *daj;
  double adaij, termj, absadajj;
  mwIndex *dzknnz, *dzir, *blksj;

  rsdpN = pcK->rsdpN;
/* ------------------------------------------------------------
   Partition working arrays
   mwIndex: dzknnz(psdN=length(K.s)), blksj(psdN), dzir(dznnz = dzstructjc[m]),
     iwork[iwsiz],
     with iwsiz = max(dznnz, max(nk(PSD))).
   double:    daj(lenud), fwork[fwsiz],
     with fwsiz = 2 * max(rMaxn^2, 2*hMaxn^2).
   ------------------------------------------------------------ */
  nblk = pcK->sdpN;
  dznnz = dzstructjc[m];
  if(dznnz <= 0)
    return;                                 /* nothing to do if no PSD*/
  dzknnz = iwork;                           /* dzknnz(nblk) */
  blksj = dzknnz + nblk;                    /* blksj(nblk) */
  dzir = blksj + nblk;                     /* dzir(dznnz) */
  iwork  = dzir + dznnz;                    /* iwork(iwsiz) */
  iwsiz -= 2*nblk + dznnz;
  mxAssert(iwsiz >= MAX(pcK->rMaxn,pcK->hMaxn), "iwork too small in getada3()");
  daj    = fwork;                           /* lenfull */
  fwork = daj + lenud;                      /* fwsiz */
  fwsiz -= lenud;
  mxAssert(fwsiz >= 2 * MAX(SQR(pcK->rMaxn),2 * SQR(pcK->hMaxn)), "fwork too small in getada3()");
/* ------------------------------------------------------------
   Make "lenfull" vector index valid into daj.
   ------------------------------------------------------------ */
  daj -= blkstart[0];     /* We'll use daj[blkstart[0]:end] */
/* ------------------------------------------------------------
   Initialize dzknnz = 0, meaning dz=[]. Later we will merge
   columns from dzstruct, with dz, and partition into selected blocks.
   ------------------------------------------------------------ */
  mxAssert(dznnz > 0, "");  /* we know that there exist nonempty PSD: */
  for(i = 0; i < nblk; i++)
    dzknnz[i] = 0;
  for(j = 1; dzstructjc[j] == 0; j++); /* 1st nonzero PSD constraint */
/* ============================================================
   MAIN getada LOOP: loop over nodes perm(0:m-1)
   ============================================================ */
  for(--j; j < m; j++){
    permj = perm[j];
/* ------------------------------------------------------------
   Make dzir: the PSD-nonzero locations, with pointers
   to the selected PSD blocks. nz-locs = merge(dzir,dzstruct(:,j)).
   ------------------------------------------------------------ */
    i = dzstructjc[j];
    while( i < dzstructjc[j+1]){
      k = xblk[dzstructir[i]];
      knz = i;                                /* add dzstructir(knz:i-1) */
      intbsearch(&i, dzstructir, dzstructjc[j+1], blkstart[k+1]);
      mxAssert(i > knz,"");
      exmerge(dzir+dzjc[k], dzstructir+knz, dzknnz[k],i-knz,iwsiz,
              cwork,iwork);
      dzknnz[k] += i-knz;                     /* number added */
    }
/* ------------------------------------------------------------
   Compute daj = P(d)*aj = vec(D*Aj*D).
   ------------------------------------------------------------ */
    nnzbj = spsqrscale(daj,blksj,dzjc,dzir,dzknnz, udsqr,
                       At.ir,At.pr,Ajc1[permj],At.jc[permj+1],
                       blkstart, xblk, psdNL, rsdpN, fwork, iwork);
    /* iwork(max(K.s)), fwork(2 * (rMaxn^2 + 2*hMaxn^2))*/
    mxAssert(nnzbj <= nblk, ""); /* number of nz-matrix-blocks */
/* ------------------------------------------------------------
   For all i with invpermi < j:
   ada_ij = a_i'*daj.
   ------------------------------------------------------------ */
    for(inz = ada.jc[permj]; inz < ada.jc[permj+1]; inz++){
      i = ada.ir[inz];
      if(invperm[i] <= j){
        adaij = ada.pr[inz];
        if(invperm[i] < j)
          for(knz = Ajc1[i]; knz < At.jc[i+1]; knz++)
            adaij +=  At.pr[knz] * daj[At.ir[knz]];
        else{         /* diag entry: absd[j] = sum(abs(aj.*daj)) */
          absadajj = adaij;
          for(knz = Ajc1[i]; knz < At.jc[i+1]; knz++){
            termj = At.pr[knz] * daj[At.ir[knz]];
            adaij +=  termj;
            absadajj += fabs(termj);
          }
          absd[permj] = absadajj;
        }
        ada.pr[inz] = adaij;
      }
    }
/* ------------------------------------------------------------
   Set daj = all-0
   ------------------------------------------------------------ */
    for(knz = 0; knz < nnzbj; knz++){
      i = blksj[knz];
      spzeros(daj,dzir+dzjc[i],dzknnz[i]);
    }
  } /* j = 0:m-1 */
  mxAssert(dzjc[nblk-1]+dzknnz[nblk-1] == dznnz,"");
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
  mxArray *myplhs[NPAROUT];
  coneK cK;
  const mxArray *MY_FIELD;
  mwIndex lenfull, lenud, m, i, j, k, fwsiz, iwsiz, dznnz, maxadd;
  const double *permPr, *Ajc1Pr, *blkstartPr, *udsqr;
  const mwIndex *dzstructjc, *dzstructir;
  double *fwork, *absd;
  mwIndex *blkstart, *iwork, *Ajc1, *psdNL, *xblk, *perm, *invperm, *dzjc;
  bool *cwork;
  jcir At, ada;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "getADA requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "getADA produces less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Compute some statistics based on cone K structure
   ------------------------------------------------------------ */
  lenud = cK.rDim + cK.hDim;                  /* for PSD */
  lenfull = cK.lpN +  cK.qDim + lenud;
/* ------------------------------------------------------------
   Allocate working array blkstart(|K.s|+1).
   ------------------------------------------------------------ */
  blkstart = (mwIndex *) mxCalloc(cK.sdpN + 1, sizeof(mwIndex));
/* ------------------------------------------------------------
   Translate blkstart from Fortran-double to C-mwIndex
   ------------------------------------------------------------ */
  MY_FIELD = mxGetField(K_IN,(mwIndex)0,"blkstart");        /*K.blkstart*/
  mxAssert( MY_FIELD != NULL, "Missing K.blkstart.");
  mxAssert(mxGetM(MY_FIELD) * mxGetN(MY_FIELD) == 2+cK.lorN+cK.sdpN, "Size mismatch K.blkstart.");
  blkstartPr = mxGetPr(MY_FIELD) + cK.lorN + 1;          /* point to start of PSD */
  for(i = 0; i <= cK.sdpN; i++){                            /* to integers */
    j = (mwIndex) blkstartPr[i];
    mxAssert(j>0,"");
    blkstart[i] = --j;
  }
/* ------------------------------------------------------------
   INPUT sparse constraint matrix At:
   ------------------------------------------------------------ */
  mxAssert(mxGetM(AT_IN) == lenfull, "Size mismatch At");          /* At */
  m = mxGetN(AT_IN);
  mxAssert(mxIsSparse(AT_IN), "At should be sparse.");
  At.pr = mxGetPr(AT_IN);
  At.jc = mxGetJc(AT_IN);
  At.ir = mxGetIr(AT_IN);
/* ------------------------------------------------------------
   Get SCALING VECTOR: udsqr
   ------------------------------------------------------------ */
  mxAssert(mxGetM(UDSQR_IN) * mxGetN(UDSQR_IN) == lenud, "udsqr size mismatch.");   /* udsqr */
  udsqr = mxGetPr(UDSQR_IN);
/* ------------------------------------------------------------
   Get Ajc1
   ------------------------------------------------------------ */
  mxAssert(mxGetM(AJC1_IN)*mxGetN(AJC1_IN) == m, "Ajc1 size mismatch");
  Ajc1Pr = mxGetPr(AJC1_IN);
/* ------------------------------------------------------------
   DISASSEMBLE Aord structure: Aord.{dz,sperm}
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(AORD_IN), "Aord should be a structure.");
  MY_FIELD = mxGetField(AORD_IN,(mwIndex)0,"dz");         /* Aord.dz */
  mxAssert( MY_FIELD != NULL, "Missing field Aord.dz.");
  mxAssert(mxGetN(MY_FIELD) >= m, "Size mismatch Aord.dz.");
  mxAssert(mxGetM(MY_FIELD) == lenfull, "Aord.dz size mismatch");
  mxAssert(mxIsSparse(MY_FIELD), "Aord.dz should be sparse.");
  dzstructjc = mxGetJc(MY_FIELD);
  dzstructir = mxGetIr(MY_FIELD);
  MY_FIELD = mxGetField(AORD_IN,(mwIndex)0,"sperm");   /* Aord.sperm */
  mxAssert( MY_FIELD != NULL, "Missing field Aord.sperm.");
  mxAssert(mxGetM(MY_FIELD) * mxGetN(MY_FIELD) == m, "Aord.sperm size mismatch");
  permPr = mxGetPr(MY_FIELD);
/* ------------------------------------------------------------
   Allocate output matrix ADA as a duplicate of ADA_IN:
   ------------------------------------------------------------ */
  mxAssert(mxGetM(ADA_IN) == m && mxGetN(ADA_IN) == m, "Size mismatch ADA.");
  mxAssert(mxIsSparse(ADA_IN), "ADA should be sparse.");
  ADA_OUT = mxDuplicateArray(ADA_IN);                   /* ADA = ADA_IN */
  ada.jc = mxGetJc(ADA_OUT);
  ada.ir = mxGetIr(ADA_OUT);
  ada.pr = mxGetPr(ADA_OUT);
/* ------------------------------------------------------------
   Create output vector absd(m)
   ------------------------------------------------------------ */
  ABSD_OUT = mxCreateDoubleMatrix(m,(mwSize)1,mxREAL);
  absd = mxGetPr(ABSD_OUT);
/* ------------------------------------------------------------
   The following ONLY if there are PSD blocks:
   ------------------------------------------------------------ */
  if(cK.sdpN > 0){
    maxadd = dzstructjc[1];
    for(i = 1; i < m; i++)
      if(dzstructjc[i+1] > dzstructjc[i] + maxadd)
        maxadd = dzstructjc[i+1] - dzstructjc[i];
/* ------------------------------------------------------------
   ALLOCATE integer work array iwork(iwsiz), with
   iwsiz = MAX(m, 2*nblk + dznnz +
     max(maxadd+2+log_2(1+maxadd), max(nk(PSD)))),
   where dznnz = dzstructjc[m].
   ------------------------------------------------------------ */
    dznnz = dzstructjc[m];
    iwsiz = (mwIndex) floor(log(1.0+maxadd)/log(2.0));  /* double to mwIndex */
    iwsiz += maxadd + 2;
    iwsiz = 2*cK.sdpN + dznnz + MAX(iwsiz,MAX(cK.rMaxn,cK.hMaxn));
    iwork = (mwIndex *) mxCalloc(MAX(iwsiz,m), sizeof(mwIndex));
/* ------------------------------------------------------------
   ALLOCATE integer working arrays:
   Ajc1(m) psdNL[cK.sdpN], dzjc(cK.sdpN+1), perm(m), invperm(m), xblk(lenud).
   cwork(maxadd).
   ------------------------------------------------------------ */
    Ajc1 = (mwIndex *) mxCalloc(MAX(m,1), sizeof(mwIndex));
    psdNL = (mwIndex *) mxCalloc(1+2*cK.sdpN + lenud, sizeof(mwIndex));
    xblk = psdNL + cK.sdpN;    /* Not own alloc: we'll subtract blkstart[0] */
    dzjc = xblk + lenud;       /*dzjc(sdpN+1) */
    perm = (mwIndex *) mxCalloc(MAX(2 * m,1), sizeof(mwIndex));
    invperm = perm + m;                                 /* invperm(m) */
    cwork = (bool *) mxCalloc(MAX(1,maxadd), sizeof(bool));
/* ------------------------------------------------------------
   ALLOCATE float working array:
   fwork[fwsiz] with fwsiz = lenud + 2 * max(rMaxn^2, 2*hMaxn^2).
   ------------------------------------------------------------ */
    fwsiz = lenud + 2 * MAX(SQR(cK.rMaxn),2*SQR(cK.hMaxn));
    fwork  = (double *) mxCalloc(MAX(fwsiz,1), sizeof(double));
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
   Let psdNL = K.s in integer, Ajc1 = Ajc1Pr in integer.
   ------------------------------------------------------------ */
    for(i = 0; i < cK.sdpN; i++)                /* K.s */
      psdNL[i] = (mwIndex) cK.sdpNL[i];
    for(i = 0; i < m; i++)
      Ajc1[i] = (mwIndex) Ajc1Pr[i];
/* ------------------------------------------------------------
   Let k = xblk(j-blkstart[0]) iff
   blkstart[k] <= j < blkstart[k+1], k=0:nblk-1.
   ------------------------------------------------------------ */
    j = blkstart[0];
    xblk -= j;     /* Make blkstart[0]:blkstart[end] valid indices */
    for(k = 0; k < cK.sdpN; k++){
      i = blkstart[k+1];
      while(j < i)
        xblk[j++] = k;
    }
/* ------------------------------------------------------------
   ACTUAL COMPUTATION: handle constraint aj=At(:,perm(j)), j=0:m-1.
   ------------------------------------------------------------ */
    dzblkpartit(dzjc, dzstructir, xblk, dznnz, cK.sdpN);
    getada3(ada, absd, At,udsqr,Ajc1, dzjc,dzstructjc,dzstructir, blkstart,
            xblk,psdNL, perm,invperm, m,lenud, &cK, fwork,fwsiz,
            iwork,iwsiz, cwork);
/* ------------------------------------------------------------
   RELEASE WORKING ARRAYS (for PSD blocks only).
   ------------------------------------------------------------ */
    mxFree(fwork);
    mxFree(cwork);
    mxFree(perm);
    mxFree(psdNL);
    mxFree(Ajc1);
  }  /* ~isempty(K.s) */
/* ------------------------------------------------------------
   If no PSD-blocks, than we merely compute absd = diag(ADA)
   ALLOCATE integer work array iwork(m), with
   ------------------------------------------------------------ */
  else{
    iwork = (mwIndex *) mxCalloc(MAX(1,m), sizeof(mwIndex));
    cpspdiag(absd, ada,m);
  }
/* ------------------------------------------------------------
   Let ADA = (ADA+ADA')/2, so that it gets symmetric.
   ------------------------------------------------------------ */
  spmakesym(ada,m,iwork);         /* uses iwork(m) */
/* ------------------------------------------------------------
   RELEASE WORKING ARRAYS iwork and blkstart.
   ------------------------------------------------------------ */
  mxFree(iwork);
  mxFree(blkstart);
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  i = MAX(nlhs, 1);
  memcpy(plhs,myplhs, i * sizeof(mxArray *));
  for(; i < NPAROUT; i++)
    mxDestroyArray(myplhs[i]);
}
