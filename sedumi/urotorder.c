/*
%                         [u,perm,gjc,g] = urotorder(u,K, maxu,permIN)
% UROTORDER  Stable reORDERing of triu U-factor by Givens ROTations.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function [u,perm,gjc,g] = urotorder(u,K, maxu,permIN)

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
#include <math.h>
#include "mex.h"
#include "blksdp.h"
#include "givens.h"

#define U_OUT myplhs[0]
#define PERM_OUT myplhs[1]
#define GJC_OUT myplhs[2]
#define G_OUT myplhs[3]
#define NPAROUT 4

#define U_IN prhs[0]
#define K_IN prhs[1]
#define MAXU_IN prhs[2]
#define NPARINMIN 3
#define PERM_IN prhs[3]
#define NPARIN 4

/* ============================================================
   TYPE DEFINITIONS:
   ============================================================ */
/* controls when d's are recomputed from scratch. Not critical, since
   d's are used only for SELECTING pivots. Just need avoid full underflow. */
#define DRELTOL 1E-10

/* ************************************************************
   PROCEDURE rotorder
   UPDATED
     u - full n x n matrix. On input, triu(u) is possibly unstable factor.
        On output, triu(u(:,perm)) is a stable factor. U_OUT = Q*U_IN,
        where Q is a sequence of givens rotations, given in g.
   OUTPUT
     perm - length n stable (column) pivot ordering.
     gjc - The givens rotations at step k are g[gjc[k]:gjc[k+1]-1].
       The order in each column is bottom up.
     g - length gjc[n] <= n(n-1)/2 array of givens rotations.
      At worst we need n-1-k rotations in iter k=0:n-2.
   ************************************************************ */
void rotorder(mwIndex *perm, double *u, mwIndex *gjc, twodouble *g, double *d,
              const double maxusqr, const mwIndex n)
{
  mwIndex i,j,k,inz, pivk, m;
  double *uj, *rowuk;
  double dk,y,nexty, h, uki,ukmax;
  twodouble gi;
/* ------------------------------------------------------------
   Initialize:
   Let perm = 1:n, inz = 0. (inz points into rotation list r)
   Let d(0) = 0, h = 1: this will let us compute all d's (since d(0)<h).
   ------------------------------------------------------------ */
  for(j = 0; j < n; j++)
    perm[j] = j;
  inz = 0;
  d[0] = 0.0; h = 1.0;
  for(k = 0, rowuk = u; k < n-1; k++, rowuk++){
    gjc[k] = inz;
/* ------------------------------------------------------------
   If current d's are not reliable then
   compute d(i) = sum(u(k:n-1,i).^2) from scratch.
   ------------------------------------------------------------ */
    if(d[perm[k]] <= h){
      for(j = k; j < n; j++){
        i = perm[j];
        d[i] = realssqr(rowuk + i*n,j+1-k);
      }
      h = d[perm[k]] * DRELTOL;
    }
/* ------------------------------------------------------------
   Let ukmax = max(U(k,perm(k+1:n)).^2)
   ------------------------------------------------------------ */
    ukmax = 0.0;
    for(j = k + 1; j < n; j++){
      uki = rowuk[perm[j] * n];
      uki *= uki;
      ukmax = MAX(ukmax, uki);
    }
/* ------------------------------------------------------------
   If ukmax >  maxusqr * d(k), then pivot k is unstable.
   If so, find best pivot: (pivk, dk) = max(perm(d(k:n))).
   ------------------------------------------------------------ */
    if(ukmax > maxusqr * d[perm[k]]){
      dk = 0.0;
      for(j = k+1; j < n; j++)
        if(d[perm[j]] > dk){
          pivk = j;
          dk = d[perm[j]];
        }
/* ------------------------------------------------------------
   Pivot on column pivk, and make U(:,perm)
   upper-triangular by pivk - k givens rotations on U(:,perm(k:n)).
   Givens at row i is {u(i,j), norm( u(i+1:pivk,j) )} for
   j=perm[pivk] and i = k:pivk-1.
   ------------------------------------------------------------ */
      m = pivk - k;                    /* number of Givens rotations needed */
      j = perm[pivk];                  /* uj(1:m) should become 0 */
      uj = rowuk + j * n;
      nexty = uj[m];                   /* last nonzero in col uj */
      y = SQR(nexty);
      for(i = m; i > 0; i--){
        gi.x = uj[i-1];
        gi.y = nexty;
        y += SQR(gi.x);
        nexty = sqrt(y);
        gi.x /= nexty;                  /* Normalize to rotation [x,y; y,-x] */
        gi.y /= nexty;
        g[i-1] = gi;
      }                                /* y == d[j] after loop */
      uj[0] = nexty;                   /* New pivotal diagonal entry */
/* ------------------------------------------------------------
   move pivot j=perm[pivk] to head of perm (shifting old k:pivk-1)
   ------------------------------------------------------------ */
      memmove(perm+k+1, perm+k, m * sizeof(mwIndex));     /* move 1-> */
      perm[k] = j;                     /* inserted at k */
/* ------------------------------------------------------------
   Apply rotations to columns perm(k+1:n-1).
   Apply 1,2,...,m rotations on column k+1,..,k+m=pivk,
   and m rotations on cols pivk+1:n-1.
   ------------------------------------------------------------ */
      for(i = 1; i <= m; i++)
        givensrotuj(rowuk + perm[k+i] * n, g,i);
      for(i += k; i < n; i++)
        givensrot(rowuk + perm[i] * n, g,m);
      inz += m;                         /* point to next avl. place */
      g += m;
/* ------------------------------------------------------------
   Update d(perm(k+1:n)) -= u(k,perm(k+1:n)).^2.
   ------------------------------------------------------------ */
      for(j = k + 1; j < n; j++){
        i = perm[j];
        d[i] -= SQR(rowuk[i * n]);
      }
    }
  }
/* ------------------------------------------------------------
   We have reordered n-1 columns of U using inz Givens-rotations.
   ------------------------------------------------------------ */
  mxAssert(n > 0,"");
  gjc[n-1] = inz;
}

/* ************************************************************
   PROCEDURE prpirotorder
   UPDATED
     u,upi - full n x n matrix. On input, triu(u) is possibly unstable factor.
        u is triu, real diagonal.
        On output, triu(u(:,perm)) is a stable factor. U_OUT = Q*U_IN,
        where Q is a sequence of givens rotations, given in g.
        u remains triu, real diagonal.
   OUTPUT
     perm - length n stable (column) pivot ordering.
     gjc - The givens rotations at step k are g[gjc[k]:gjc[k+1]-1].
       The order in each column is bottom up.
     g - length gjc[n] <= n(n-1)/2 array of givens rotations.
      At worst we need n-1-k rotations in iter k=0:n-2.
   ************************************************************ */
void prpirotorder(mwIndex *perm, double *u,double *upi, mwIndex *gjc,
                  tridouble *g, double *d,
                  const double maxusqr, const mwIndex n)
{
  mwIndex i,j,k,inz, pivk, m;
  double *uj,*ujpi, *rowuk, *rowukpi;
  double dk,y,nexty, h, uki,ukiim,ukmax;
  tridouble gi;
/* ------------------------------------------------------------
   Initialize:
   Let perm = 1:n, inz = 0. (inz points into rotation list r)
   Let d(0) = 0, h = 1: this will let us compute all d's (since d(0)<h).
   ------------------------------------------------------------ */
  for(j = 0; j < n; j++)
    perm[j] = j;
  inz = 0;
  d[0] = 0.0; h = 1.0;
  for(k = 0, rowuk = u, rowukpi = upi; k < n-1; k++, rowuk++, rowukpi++){
    gjc[k] = inz;
/* ------------------------------------------------------------
   If current d's are not reliable then
   compute d(i) = sum(u(k:n-1,i).^2) from scratch.
   ------------------------------------------------------------ */
    if(d[perm[k]] <= h){
      for(j = k; j < n; j++){
        i = perm[j];                /* diag entry u(j,i) is real */
        d[i] = realssqr(rowuk + i*n,j+1-k) + realssqr(rowukpi + i*n,j-k);
      }
      h = d[perm[k]] * DRELTOL;
    }
/* ------------------------------------------------------------
   Let ukmax = max(abs(U(k,perm(k+1:n))).^2)
   ------------------------------------------------------------ */
    ukmax = 0.0;
    for(j = k + 1; j < n; j++){
      uki = rowuk[perm[j] * n];
      ukiim = rowukpi[perm[j] * n];
      ukmax = MAX(ukmax, SQR(uki) + SQR(ukiim));
    }
/* ------------------------------------------------------------
   If ukmax >  maxusqr * d(k), then pivot k is unstable.
   If so, find best pivot: (pivk, dk) = max(perm(d(k:n))).
   ------------------------------------------------------------ */
    if(ukmax > maxusqr * d[perm[k]]){
      dk = 0.0;
      for(j = k+1; j < n; j++)
        if(d[perm[j]] > dk){
          pivk = j;
          dk = d[perm[j]];
        }
/* ------------------------------------------------------------
   Pivot on column pivk, and make U(:,perm)
   upper-triangular by pivk - k givens rotations on U(:,perm(k:n)).
   Givens at row i is {u(i,j), norm( u(i+1:pivk,j) )} for
   j=perm[pivk] and i = k:pivk-1.
   Thus, the rotation gi consists of 1 complex and 1 real entry,
   (gi.x,gi.xim) and gi.y, resp. The rotation is [conj(x), y;y, -x]
   ------------------------------------------------------------ */
      m = pivk - k;                    /* number of Givens rotations needed */
      j = perm[pivk];                  /* uj(1:m) should become 0 */
      uj = rowuk + j * n;
      ujpi = rowukpi + j * n;
      nexty = uj[m];                   /* last nonzero in col uj (real) */
      y = SQR(nexty);
      for(i = m; i > 0; i--){
        gi.x = uj[i-1];
        gi.xim = ujpi[i-1];
        gi.y = nexty;
        y += SQR(gi.x) + SQR(gi.xim);
        nexty = sqrt(y);
        gi.x /= nexty;         /* Normalize to rotation [conj(x),y; y,-x] */
        gi.xim /= nexty;
        gi.y /= nexty;
        g[i-1] = gi;
      }                                /* y == d[j] after loop */
      uj[0] = nexty;                   /* New pivotal diagonal entry */
/* ------------------------------------------------------------
   move pivot j=perm[pivk] to head of perm (shifting old k:pivk-1)
   ------------------------------------------------------------ */
      memmove(perm+k+1, perm+k, m * sizeof(mwIndex));     /* move 1-> */
      perm[k] = j;                     /* inserted at k */
/* ------------------------------------------------------------
   Apply rotations to columns perm(k+1:n-1).
   Apply 1,2,...,m rotations on column k+1,..,k+m=pivk,
   and m rotations on cols pivk+1:n-1.
   ------------------------------------------------------------ */
      for(i = 1; i <= m; i++)
        prpigivensrotuj(rowuk + perm[k+i] * n,rowukpi + perm[k+i] * n, g,i);
      for(i += k; i < n; i++)
        prpigivensrot(rowuk + perm[i] * n,rowukpi + perm[i] * n, g,m);
      inz += m;                         /* point to next avl. place */
      g += m;
/* ------------------------------------------------------------
   Update d(perm(k+1:n)) -= u(k,perm(k+1:n)).^2.
   ------------------------------------------------------------ */
      for(j = k + 1; j < n; j++){
        i = perm[j];
        d[i] -= SQR(rowuk[i * n]) + SQR(rowukpi[i * n]);
      }
    }
  }
/* ------------------------------------------------------------
   We have reordered n-1 columns of U using inz Givens-rotations.
   ------------------------------------------------------------ */
  mxAssert(n > 0,"");
  gjc[n-1] = inz;
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  mxArray *myplhs[NPAROUT];
  mwIndex i,j,k, nk, nksqr, lenud, sdplen, gnnz, inz, maxKs,maxKssqr, rgnnz, hgnnz;
  const double *uOld, *permOld;
  double *u, *d, *gjcPr, *permPr, *fwork, *fworkpi;
  mwIndex *perm, *gjc;
  double *g, *gk;
  double maxusqr;
  coneK cK;
  char use_pivot;
/* ------------------------------------------------------------
   Check for proper number of arguments 
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARINMIN, "urotorder requires more input arguments.");
  mxAssert(nlhs <= NPAROUT, "urotorder generates less output arguments.");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Get statistics of cone K structure
   ------------------------------------------------------------ */
  lenud = cK.rDim + cK.hDim;
  sdplen = cK.rLen + cK.hLen;
/* ------------------------------------------------------------
   Get scalar input MAXU and input vectors U_IN, PERM_IN
   ------------------------------------------------------------ */
  maxusqr = mxGetScalar(MAXU_IN);
  maxusqr *= maxusqr;
  mxAssert(mxGetM(U_IN) * mxGetN(U_IN) == lenud, "u size mismatch");
  uOld = mxGetPr(U_IN);
  use_pivot = 0;
  if(nrhs >= NPARIN)                              /* Optional permIN */
    if(mxGetM(PERM_IN) * mxGetN(PERM_IN) > 0){
      mxAssert(mxGetM(PERM_IN) * mxGetN(PERM_IN) == sdplen, "perm size mismatch");
      use_pivot = 1;
      permOld = mxGetPr(PERM_IN);
      }
/* ------------------------------------------------------------
   Allocate output U_OUT, and initialize u_out = u_in.
   ------------------------------------------------------------ */
  U_OUT = mxCreateDoubleMatrix(lenud, (mwSize)1, mxREAL);
  u = mxGetPr(U_OUT);
  memcpy(u, mxGetPr(U_IN), lenud * sizeof(double));
/* ------------------------------------------------------------
   Allocate outputs PERM(sum(K.s)), GJC(sum(K.s))
   ------------------------------------------------------------ */
  PERM_OUT = mxCreateDoubleMatrix(sdplen, (mwSize)1, mxREAL);
  permPr = mxGetPr(PERM_OUT);
  GJC_OUT  = mxCreateDoubleMatrix(sdplen, (mwSize)1, mxREAL);
  gjcPr = mxGetPr(GJC_OUT);
/* ------------------------------------------------------------
   Allocate g initially as length (lenud - cK.rLen) / 2. The final
   length can be shorter (viz. gjc[sum(K.s)])
   ------------------------------------------------------------ */
  rgnnz = (cK.rDim - cK.rLen) / 2;               /* n(n-1)/2 real sym */
  hgnnz = (cK.hDim - 2*cK.hLen) / 4;             /* n(n-1)/2 complex herm */
  gnnz = rgnnz * 2 + hgnnz * 3;
  g = (double *) mxCalloc(MAX(1, gnnz),sizeof(double));
/* ------------------------------------------------------------
   Allocate working arrays:
   Let maxKssqr = max(rMaxn^2, 2*hMaxn^2), then
   integer perm(max(K.s)), gjc(max(K.s))
   double d(max(K.s)), fwork(maxKs)
   ------------------------------------------------------------ */
  maxKs = MAX(cK.rMaxn,cK.hMaxn);                     /* max(K.s) */
  maxKssqr = MAX(SQR(cK.rMaxn),2 * SQR(cK.hMaxn));    /* max(K.s.^2) */
  perm = (mwIndex *) mxCalloc(MAX(1,maxKs), sizeof(mwIndex));
  gjc  = (mwIndex *) mxCalloc(MAX(1,maxKs), sizeof(mwIndex));
  d     = (double *) mxCalloc(MAX(1,maxKs), sizeof(double));
  fwork = (double *) mxCalloc(MAX(1,maxKssqr), sizeof(double));
  fworkpi = fwork + SQR(cK.hMaxn);
/* ------------------------------------------------------------
   The actual job is done here: U_NEW = Q(g) * U_OLD
   ------------------------------------------------------------ */
  inz = 0;
  for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
    nk = cK.sdpNL[k];
    nksqr = SQR(nk);
    memcpy(fwork, uOld, nksqr *sizeof(double));   /* k-th U-matrix */
    gk = g+inz;
    rotorder(perm, fwork, gjc, (twodouble *) gk, d, maxusqr, nk);
/* ------------------------------------------------------------
   Physically reorder the columns from fwork into u. Then Let
   tril(U) = triu(U)'
   ------------------------------------------------------------ */
    uperm(u, fwork, perm, nk);
    triu2sym(u,nk);
/* ------------------------------------------------------------
   Let perm_out = perm_in(perm)
   ------------------------------------------------------------ */
    if(use_pivot){
      for(i = 0; i < nk; i++)
        permPr[i] = permOld[perm[i]];
      permOld += nk;
    }
    else
      for(i = 0; i < nk; i++)
        permPr[i] = 1.0 + perm[i];
    for(i = 0; i < nk; i++)
      gjcPr[i] = gjc[i];               /* don't add 1 */
    inz += 2 * gjc[nk-1];     /* next PSD block. Rotation g is 2 doubles */
    gjcPr += nk;
    permPr += nk; uOld += nksqr;
    u += nksqr;
  }
/* ------------------------------------------------------------
   Complex Hermitian
   ------------------------------------------------------------ */
  for(; k < cK.sdpN; k++){                    /* complex Hermitian */
    nk = cK.sdpNL[k];
    nksqr = SQR(nk);
    memcpy(fwork, uOld, nksqr *sizeof(double));   /* k-th complex U-matrix */
    memcpy(fworkpi, uOld+nksqr, nksqr *sizeof(double));
    gk = g+inz;
    prpirotorder(perm, fwork,fworkpi, gjc, (tridouble *) gk, d, maxusqr, nk);
/* ------------------------------------------------------------
   Physically reorder the columns from fwork into u. Then Let
   tril(U) = triu(U)'
   ------------------------------------------------------------ */
    uperm(u, fwork, perm, nk);                  /* real part */
    uperm(u+nksqr, fworkpi, perm, nk);      /* imaginary part */
    triu2herm(u,u+nksqr, nk);
/* ------------------------------------------------------------
   Let perm_out = perm_in(perm)
   ------------------------------------------------------------ */
    if(use_pivot){
      for(i = 0; i < nk; i++)
        permPr[i] = permOld[perm[i]];
      permOld += nk;
    }
    else
      for(i = 0; i < nk; i++)
        permPr[i] = 1.0 + perm[i];
    for(i = 0; i < nk; i++)
      gjcPr[i] = gjc[i];               /* don't add 1 */
    inz += 3 * gjc[nk-1];     /* next PSD block. Rotation g is 3 doubles */
    gjcPr += nk;
    permPr += nk;
    nksqr += nksqr;
    uOld += nksqr; u += nksqr;
  }
/* ------------------------------------------------------------
   In total, we used inz doubles in Givens rotations.
   Reallocate (shrink) g accordingly.
   ------------------------------------------------------------ */
  mxAssert(inz <= gnnz,"");
  if(inz > 0){
    if((g = (double *) mxRealloc(g, inz * sizeof(double))) == NULL)
      mexErrMsgTxt("Memory allocation error.");
  }
  else{
    mxFree(g);
    g = (double *) NULL;
  }
/* ------------------------------------------------------------
   Assign g to a length inz output vector
   ------------------------------------------------------------ */
  G_OUT = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
  mxFree(mxGetPr(G_OUT));
  mxSetPr(G_OUT, (double *) g);
  mxSetM(G_OUT, inz);
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(fwork);
  mxFree(d);
  mxFree(gjc);
  mxFree(perm);
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  i = MAX(nlhs, 1);
  memcpy(plhs,myplhs, i * sizeof(mxArray *));
  for(; i < NPAROUT; i++)
    mxDestroyArray(myplhs[i]);
}
