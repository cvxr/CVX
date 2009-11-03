/*

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
#include "blksdp.h"
#include "mex.h"

/* ------------------------------------------------------------
   PROTOTYPES:
   ------------------------------------------------------------ */
mwIndex blkLDL(const mwIndex neqns, const mwIndex nsuper, const mwIndex *xsuper,
           const mwIndex *snode,  const mwIndex *xlindx, const mwIndex *lindx,
           double *lb,
           const mwIndex *ljc, double *lpr, double *d, const mwIndex *perm,
           const double ub, const double maxu, mwIndex *skipIr,
           mwIndex iwsiz, mwIndex *iwork, mwIndex fwsiz, double *fwork);

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- isscalarmul(x,alpha,n)
   Computes x *= alpha using BLAS.
   ************************************************************ */
void isscalarmul(double *x, const double alpha, const mwIndex n)
{
    mwIndex one=1;
    #ifdef PC
    dscal(&n,&alpha,x,&one);
    #endif
    #ifdef UNIX
    dscal_(&n,&alpha,x,&one);
    #endif
    return;
}

/* ************************************************************
   PROCEDURE maxabs - computes inf-norm using BLAS
   INPUT
     x - vector of length n
     n - length of x.
   RETURNS y = norm(x,inf).
   ************************************************************ */
double maxabs(const double *x,const mwIndex n)
{
mwIndex one=1;
#ifdef PC
    return fabs(x[idamax(&n,x,&one)]);
#endif
#ifdef UNIX
    return fabs(x[idamax_(&n,x,&one)]);
#endif
}

/* ************************************************************
   PROCEDURE cholonBlk - CHOLESKY on a dense diagonal block.
            Also updates nonzeros below this diagonal block -
            they need merely be divided by the scalar diagonals
            "lkk" afterwards.
   INPUT
     m      - number of rows (length of the first column).
     ncols  - number of columns in the supernode.(n <= m)
     lb     - Length ncols. Skip k-th pivot if drops below lb[k].
     ub     - max(diag(x)) / maxu^2. No stability check for pivots > ub.
     maxu   - Max. acceptable |lik|/lkk when lkk suffers cancelation.
     first - global column number of column 0. This is used only to insert
        the global column numbers into skipIr.
   UPDATED
     x  - On input, contains the columns of the supernode to
          be factored. On output, contains the factored columns of
          the supernode.
     skipIr - Lists skipped pivots with their global column number
        in 0:neqns-1. Active range is first:first+ncols-1.
        Skipped if d(k) suffers cancelation and max(abs(L(:,k)) > maxu.
     *pnskip - nnz in skip; *pnskip <= order N of sparse matrix.
   OUTPUT
     d - Length ncols. Diagonal in L*diag(d)*L' with diag(L)=all-1.
   ************************************************************ */
void cholonBlk(double *x, double *d, mwIndex m, const mwIndex ncols, const mwIndex first,
               const double ub, const double maxu, double *lb,
               mwIndex *skipIr, mwIndex *pnskip)
{
  mwIndex inz,i,k,n,coltail, nskip;
  double xkk, xik, ubk;
  double *xi;
/* ------------------------------------------------------------
   Initialize:
   ------------------------------------------------------------ */
  n = ncols;
  nskip = *pnskip;
  inz = 0;
  coltail = m - ncols;
  for(k = 0; k < ncols; k++, --m, --n){
/* -------------------------------------------------------
   Let xkk = L(k,k), ubk = max(|xik|) / maxu.
   ------------------------------------------------------- */
    xkk = x[inz];
    if(xkk > lb[k]){ /* now xkk > 0 */
      if(xkk < ub){
        ubk = maxabs(x+inz+1,m-1) / maxu;
        if(xkk < ubk){
/* ------------------------------------------------------------
   If we need to add on diagonal, store this in (skipIr, lb(k)).
   ------------------------------------------------------------ */
          skipIr[nskip++] = first + k;
	  lb[k] = ubk - xkk;           /* amount added on diagonal */
	  xkk = ubk;
	}
      }
/* --------------------------------------------------------------
   Set dk = xkk, lkk = 1 (for LDL').
   -------------------------------------------------------------- */
      d[k] = xkk;                   /* now d[k] > 0 MEANS NO-SKIPPING */
      x[inz] = 1.0;
      xi = x + inz + m;                 /* point to next column */
      ++inz;
/* --------------------------------------------------------------
   REGULAR JOB: correct remaining n-k cols with col k.
   x(k+1:m,k+1:n) -= x(k+1:m,k) * x(k+1:n,k)' / xkk
   x(k+1:n,k) /= xkk,
   -------------------------------------------------------------- */
      for(i = 1; i < n; i++){
        xik = x[inz] / xkk;
        subscalarmul(xi, xik, x+inz, m-i);
        x[inz++] = xik;
        xi += m-i;
      }
      inz += coltail;                 /* Let inz point to next column */
    }
/* ------------------------------------------------------------
   If skipping is enabled and this pivot is too small:
   1) don't touch L(k:end,k): allows pivot delaying if desired.
   2) List first+k in skipIr. Set dk = 0 (MEANS SKIPPING).
   -------------------------------------------------------------- */
    else{
      skipIr[nskip++] = first + k;
      d[k] = 0.0;                 /* tag "0": means column skipped in LDL'.*/
      inz += m;                   /* Don't touch nor use L(k:end,k) */
    }
  } /* k=0:ncols-1 */
/* ------------------------------------------------------------
   Return updated number of added or skipped pivots.
   ------------------------------------------------------------ */
  *pnskip = nskip;
}

/* ************************************************************
  getbwIrInv  --  Inverse of the subscript function: given a subscript,
      irInv yields the position, counted FROM THE BOTTOM of the sparse column.

  INPUT PARAMETERS -
     nnz    - LENGTH OF THE FIRST COLUMN OF THE SUPERNODE,
              INCLUDING THE DIAGONAL ENTRY.
     Lir    - Lir[0:nnz-1] ARE THE ROW INDICES OF THE NONZEROS
              OF THE FIRST COLUMN OF THE SUPERNODE.
  OUTPUT PARAMETERS - 
     irInv - On return, irInv[Lir[0:nnz-1]] = nnz:-1:1, so that
		           Lir[nnz-irInv[i]]  == i
             The position of subscript "xij" is thus
			   xjc[j+1] - irInv[i].
   ************************************************************ */
void getbwIrInv(mwIndex *irInv, const mwIndex *Lir, const mwIndex nnz)
{
  mwIndex inz,bwinz;

  bwinz = nnz;
  for(inz = 0; inz < nnz; inz++, bwinz--)
    irInv[Lir[inz]] = bwinz;               /* bwinz = nnz:-1:1 */
}

/* ************************************************************
  suboutprod  --  Computes update from a single previous column "xk" on
		a supernode "xj", using dense computations.
  INPUT
     mj, nj  -  supernode "xj" is mj x nj.  More precisely, the column
                lengths are {mj, mj-1, ..., mj-(nj-1)}.
     xkk     -  scalar, the 1st nj entries in xk are divided by this number.
     mk      -  length of xk.  WE ASSUME mk <= mj.  Only 1st mk rows in xj
                are updated.
  UPDATED
     xj  -  On return, xj -= xk*xk(0:nj-1)'/xkk
     xk  -  On return, xk(0:nj-1) /= xkk
   ************************************************************ */
void suboutprod(double *xj, mwIndex mj, const mwIndex nj, double *xk,
                const double xkk, mwIndex mk)
{
  mwIndex j;
  double xjk;

  for(j = 0; j < nj; j++){
    xjk = xk[0] / xkk;
    subscalarmul(xj, xjk, xk, mk);   /* xj -= xjk * xk */
    xk[0] = xjk;                     /* FINAL entry ljk */
    xj += mj;                    /* point to next column which is 1 shorter */
    --mj; --mk; ++xk;
  }
}

/* ************************************************************
  isminoutprod  --  Computes update from a column "xk" and stores it in "xj",
	       using dense computations. If "xkk<=0", then let xj = 0.
  INPUT
     mk, nj  -  output "xj" is mk x nj - nj*(nj-1)/2. Its column lengths are
	        {mk, mk-1, ..., mk-(nj-1)}.
     xkk     -  scalar, the 1st nj entries in xk are divided by this number.
  OUTPUT
     xj      -  On return, xj = -xk*xk(0:nj-1)'/xkk       (NOTE THE MINUS !)
                BUT: if xkk <= 0, then xj = zeros(nj*(2m-nj+1)/2,1).
  UPDATED
     xk      -  On return, xk(0:nj-1) /= xkk if xkk > 0, otherwise unchanged.
   ************************************************************ */
void isminoutprod(double *xj, const mwIndex nj, double *xk, const double xkk,
                  mwIndex mk)
{
  mwIndex j;
  double xjk;

  if(xkk > 0.0)   /* if not phase 2 node */
    for(j = 0; j < nj; j++){
      xjk = xk[0] / xkk;
      memcpy(xj,xk,mk * sizeof(double));
      isscalarmul(xj, -xjk, mk);          /* xj = -xjk * xk */
      xk[0] = xjk;                     /* FINAL entry ljk */
      xj += mk;                /* point to next column which is 1 shorter */
      --mk; ++xk;
    }
  else  /* initialize to all-0 if phase-2 node */
    fzeros(xj,(nj * (mk + mk-nj + 1))/2);
}

/* ************************************************************
   spsuboutprod  --  Computes update from a single previous column "xk" on
		     a supernode "xj", with a different sparsity structure.
                     The relevant nonzeros of xj are accessed by a single
                     indirection, via "relind[:]".
   INPUT
     mj, nj  -  supernode "xj" has mj rows in its 1st column. In total, we
	        will update nj columns, corresponding to the 1st nj nonzeros
                in xk.
     xjjc    - xjjc[0] is start of 1st column of xj (as index into xnz), etc.
     xkk     -  scalar, the 1st nj entries in xk are divided by this number.
     mk      -  length of xk.  WE ASSUME mk <= mj.
     relind  - (mj - relind[0:mk-1]) yields the locations in xj on which the
	       xk nonzeros will act.
  UPDATED
     xnz  -  On return, xj(relind,:) -= xk*xk(0:nj-1)'/xkk
     xk   -  On return, xk(0:nj-1) /= xkk
   ************************************************************ */
void spsuboutprod(const mwIndex *xjjc, double *xnz, const mwIndex mj, const mwIndex nj,
                  double *xk,const double xkk,const mwIndex mk, const mwIndex *relind)
{
  mwIndex i, j, jcol, bottomj;
  double xjk;

  ++xjjc;             /* now it points beyond bottom of columns */
  for(j = 0; j < nj; j++){
    jcol = mj - relind[j];       /* affected column */
    bottomj = xjjc[jcol];
    xjk = xk[j] / xkk;
    for(i = j; i < mk; i++)
      xnz[bottomj - relind[i]] -= xjk * xk[i];
    xk[j] = xjk;                     /* FINAL entry ljk */
  }
}

/* ************************************************************
  spadd  --  Let xj += xk, where the supernode "xj", has a sparsity
        structure. The relevant nonzeros of xj are accessed by a indirection,
	via "relind[:]".
  INPUT
     mj, nj  -  supernode "xj" has mj rows in its 1st column. In total, we
	        will update nj columns, corresponding to the 1st nj nonzero
                rows in xk.
     xjjc    - xjjc[0] is start of 1st column of xj (as index into xnz), etc.
     mk      -  length of xk.  WE ASSUME mk <= mj.
     relind  - (mj - relind[0:mk-1]) yields the locations in xj on which the
	       xk nonzeros will act.
     xk      -  mk * nk - nk*(nk-1)/2 matrix, with column lengths
	        mk, mk-1, mk-2,.. mk-(nj-1).  They'll be substracted from
                the entries in xj that are listed by relind.
  UPDATED
     xnz     -  On return, xj(relind,:) += xk
   ************************************************************ */
void spadd(const mwIndex *xjjc, double *xnz, const mwIndex mj, const mwIndex nj,
           const double *xk, const mwIndex mk, const mwIndex *relind)
{
  mwIndex i, j, jcol, bottomj,mkcol;

  ++xjjc;             /* now it points beyond bottom of columns */
  mkcol = mk;         /* mkcol = mk - j */
  for(j = 0; j < nj; j++){
    jcol = mj - relind[j];       /* affected column */
    bottomj = xjjc[jcol];
    for(i = j; i < mk; i++)
      xnz[bottomj - relind[i]] += xk[i];
    xk += (--mkcol);   /* xk(i:mk-1) is next column */
  }
}

/* ************************************************************      
   PROCEDURE precorrect  -  Apply corrections from affecting supernode
      (skipping subnodes with non-positive diagonal) on supernodal
      diagonal block in L-factor.
   INPUT
     ljc   - start of columns in lpr.
     d     - Length neqns vector. The diagonal in L*diag(d)*L'. Only
       d[firstk:nextk-1] will be used.
     irInv - For row-indices Jir of affected supernode, Jir[m-irInv[i]]  == i.
     nextj - Last subnode of affected supernode is nextj-1.
     firstk, nextk - subnodes of affecting supernode are firstk:nextk-1.
     Kir   - unfinished row indices of affecting supernode
     mk    - number of unfinished nonzeros in affecting supernode
     fwsiz - Allocated length of fwork.
   UPDATED
     lpr  - For each column k=firstk:nextk-1, and the affected columns j
       in node, DO  L(:,j) -= (ljk / lkk) * L(:,k),
       and store the definitive j-th row of L, viz. ljk /= lkk.
   WORKING ARRAYS
     relind - length mk integer array
     fwork  - length fwsiz vector, for storing -Xk * inv(LABK) * Xk'.
   RETURNS  ncolup, number of columns updated by snode k.
    if -1, then fwsiz is too small.
   ************************************************************ */
mwIndex precorrect(double *lpr, const mwIndex *ljc,const double *d, const mwIndex *irInv,
               const mwIndex nextj, const mwIndex *Kir, const mwIndex mk,
               const mwIndex firstk, const mwIndex nextk,
               mwIndex *relind, const mwIndex fwsiz, double *fwork)
{
  mwIndex i,j,k,ncolup,mj;
  double *xj;
/* ------------------------------------------------------------
   j = first subscript in k (== 1st affected column)
   i = last subscript in k
   ncolup = number of nz-rows in k corresponding to columns in node.
   mj = number of nonzeros in l(:,j), the 1st affected column
   ------------------------------------------------------------ */
  j = Kir[0];
  i = Kir[mk-1];
  if(i < nextj)
    ncolup = mk;
  else
    for(ncolup = 1; Kir[ncolup] < nextj; ncolup++);
  mj = ljc[j+1] - ljc[j];
/* ------------------------------------------------------------
   If nz-structure of k is a single block in structure of node,
   (i.e. irInv[Kir[0]] - irInv[Kir[mk-1]] == mk-1). The subnodes
   of "node" must then be consecutive and at the start.
   Thus, we use dense computations :
   ------------------------------------------------------------ */
  if(irInv[j] - irInv[i] < mk){
    xj = lpr + ljc[j];
    for(k = firstk; k < nextk; k++)
      if(d[k] > 0.0)                        /* Skip pivot when d[k] <= 0 */
        suboutprod(xj, mj, ncolup, lpr + ljc[k+1] - mk, d[k], mk);
  }
  else{
/* ------------------------------------------------------------
   Otherwise, the nz-indices of k are scattered within the structure of node.
   Let relind be the position of these nz's in node, COUNTED FROM THE BOTTOM.
   ------------------------------------------------------------*/
    for(i = 0; i < mk; i++)
      relind[i] = irInv[Kir[i]];
/* ------------------------------------------------------------
   If k is a single column, then perform update directly in lpr:
   ------------------------------------------------------------ */
    if(nextk - firstk == 1){
      if(d[firstk] > 0.0)                   /* Skip pivot when d[k] <= 0 */
        spsuboutprod(ljc+j,lpr,mj, ncolup, lpr + ljc[nextk]-mk,
                     d[firstk],mk, relind);
    }
    else{
/* ------------------------------------------------------------
   Multiple columns in affecting snode:
   1. compute the complete modification, and store it in fwork:
   fwork = -Xk * inv(LABK) * Xk'
   ------------------------------------------------------------ */
      if(fwsiz + ncolup*(ncolup-1)/2 < mk * ncolup )
        return (mwIndex)-1;
      for(k = firstk; k < nextk; k++)      /* find 1st positive diag */
        if(d[k] > 0.0)
          break;
      if(k < nextk){                       /* if any positive diag: */
        isminoutprod(fwork, ncolup, lpr + ljc[k+1] - mk, d[k], mk);
        for(++k; k < nextk; k++)           /* remaining cols */
          if(d[k] > 0.0)                   /* Skip pivot when d[k] <= 0 */
            suboutprod(fwork, mk, ncolup, lpr + ljc[k+1] - mk, d[k], mk);
/* ------------------------------------------------------------
   2. subtract fwork from the sparse columns of node, using relind.
   ------------------------------------------------------------ */
        spadd(ljc+j,lpr,mj, ncolup, fwork,mk, relind);
      } /* end exists positive diag */
    } /* end multiple affecting cols */
  } /* end of scattered case */
/* ------------------------------------------------------------
   RETURN number of columns updated, i.e. #subnodes in k that we finished.
   ------------------------------------------------------------ */
  return ncolup;
}

/* ************************************************************
   BLKLDL  --  Block-sparse L*D*L' Cholesky factorization.

   INPUT:
      neqns   - Order "m": L is neqns * neqns
      nsuper  - Number of supernodes (blocks).
      xsuper  - Length nsuper+1: first simple-node of each supernode
      snode   - Length neqns: snode(node) is the supernode containing "node".
      xlindx  - Length nsuper+1: Start of sparsity structure in lindx,
              for each supernode (all simple nodes in a supernode have the
              same nonzero-structure).
      lindx   - row indices, for each supernode.
      ljc     - Length neqns+1: start of the columns of L.
      perm    - Length neqns: reordering of pne->At columns in Cholesky.
      ub     - max(diag(x)) / maxu^2. No stability check for pivots > ub.
      maxu    - Force max(max(abs(L))) <= maxu (by adding low-rank diag).
      iwsiz, fwsiz - size of integer and floating-point working storage.
               See "WORKING ARRAYS" for required amount.
   UPDATED:
      Lpr     - On input, contains tril(X), on output, L is
	      such that   X = L*D*L'. For columns k where d[k]=0, L(:,k)
              contains the column updated upto pivot k-1.
      lb    - Length neqns. INPUT: cancelation threshold per pivot. Skip pivot
          if it drops below.
          OUTPUT: lb(skipIr) are values of low rank diag. matrix that is
          added before factorization.
   OUTPUT
      d      - length neqns vector, diagonal in L*diag(d)*L'.
      skipIr - length nskip (<= neqns) array. skipIr(1:nskip) lists the
        columns that have been skipped in the Cholesky. d[skipIr] = 0.
   WORKING ARRAYS:
      iwork  - Length iwsiz working array, used for
           link(nsuper), length(nsuper),
           irInv(neqns), relind(neqns),
           iwsiz = 2*m + 2 * nsuper
      fwork  - Length fwsiz. Used for fwork(L.tmpsiz) in precorrect.
           fwsiz = L.tmpsiz.
   ACKNOWLEDGMENT:
       Parts are inspired by F77-block Cholesky of Ng and Peyton (ORNL).
   RETURNS  nskip (<=neqns), number of skipped nodes. Length of skipIr.
     if -1 then not enough workspace (iwsiz, fwsiz) allocated.
   ************************************************************ */
mwIndex blkLDL(const mwIndex neqns, const mwIndex nsuper, const mwIndex *xsuper,
           const mwIndex *snode,  const mwIndex *xlindx, const mwIndex *lindx,
           double *lb,
           const mwIndex *ljc, double *lpr, double *d, const mwIndex *perm,
           const double ub, const double maxu, mwIndex *skipIr,
           mwIndex iwsiz, mwIndex *iwork, mwIndex fwsiz, double *fwork)
{
  const mwIndex *Jir;
  mwIndex *link, *length, *irInv, *relind, *ncolupLst;
  mwIndex node,nextj,i,j,nnzj,n,  k,mk,linkk, snodei, nskip;
/* ------------------------------------------------------------
   Partition integer working array of size 2*(nsuper+neqns):
   iwork = [link(nsuper); length(nsuper); irInv(neqns); relind(neqns)].
   ------------------------------------------------------------ */
  if(iwsiz < 2 * (neqns + nsuper))
    return (mwIndex)-1;
  link   = iwork;                    /* 2 times length nsuper: */
  length = link + nsuper;
  irInv  = length + nsuper;          /* 2 * length neqns: */
  relind    = irInv + neqns;
/* ------------------------------------------------------------
   ncolupLst(neqns) shares the same working array as irInv(neqns).
   Namely, at stage j=xsuper[node], irInv uses only entries >= j,
   whereas ncolupLst only applies to the "old" columns < j.
   ------------------------------------------------------------ */
  ncolupLst = irInv;
/* ------------------------------------------------------------
   Initialize: link = nsuper * ones(nsuper,1) (means END-OF-LIST)
   ------------------------------------------------------------ */
  for(node = 0; node < nsuper; node++)
    link[node] = nsuper;
/* ------------------------------------------------------------
   Initialize nskip = 0.
   ------------------------------------------------------------ */
  nskip = 0;
/* ------------------------------------------------------------
   For each supernode "node", start at subnode j = xsuper[node],
   having sparsity pattern Jir.
   ------------------------------------------------------------ */
  nextj = xsuper[0];
  for(node = 0; node < nsuper; node++){
    j = nextj;		                /* 1st col in node */
    nextj = xsuper[node+1];
    n = nextj - j;			/* length of node */
    Jir = lindx + xlindx[node];         /* row-indices for column j */
    nnzj = ljc[j+1] - ljc[j];           /* nnz( column j ) */
/* ------------------------------------------------------------
   Compute inverse of Jir, yielding position from the bottom:
   Jir[nnzj-irInv[i]]  == i
   This will be handy when adding a column with a sub-sparsity structure
   to column j.
   ------------------------------------------------------------ */
    getbwIrInv(irInv, Jir, nnzj);
/* ------------------------------------------------------------
   Apply corrections from relevant previous super-nodes;
   these snodes are
   node -> link[node] -> link[link[node]] -> ...
   ------------------------------------------------------------ */
    for(k = link[node]; k < nsuper; k = link[k]){
      if((ncolupLst[k] = precorrect(lpr,ljc,d,irInv, nextj,
                                    lindx + xlindx[k+1]-length[k],
                                    length[k],xsuper[k],xsuper[k+1],
                                    relind,fwsiz,fwork)) == (mwIndex)-1 )
        return (mwIndex)-1;         /* fwsiz too small */
    }	
/* ------------------------------------------------------------
   DO DENSE CHOLESKY on the current supernode
   ------------------------------------------------------------ */
    cholonBlk(lpr + ljc[j],d+j, nnzj, n, j, ub, maxu, lb+j, skipIr,&nskip);
/* ------------------------------------------------------------
   insert each current affecting snode k into linked list of
   next supernode it will affect.
   ------------------------------------------------------------ */
    for(k = link[node]; k < nsuper; k = linkk){
      linkk = link[k];
      mk = (length[k] -= ncolupLst[k]);    /* unfinished nonzeros in k */
      if(mk){                              /* if not yet terminated: */
        i = lindx[xlindx[k+1]-mk];
        snodei = snode[i];
        link[k] = link[snodei];            /* prev. also affecting i */
        link[snodei] = k;                  /* next snode it'll affect */
      }
    }
/* ------------------------------------------------------------
   The same for current snode "node" itself:
   ------------------------------------------------------------ */
    if((length[node] = nnzj - n) > 0){
      i = Jir[n];                    /* 1st row outside snode */
      snodei = snode[i];
      link[node] = link[snodei];     /* prev. also affecting i */
      link[snodei] = node;
    }
    else
      length[node] = 0;              /* Supernode terminated */
  } /* node = 0:nsuper-1 */
/* ------------------------------------------------------------
   FINISHING: return the number of skipped pivots
   ------------------------------------------------------------ */
  return nskip;
}
