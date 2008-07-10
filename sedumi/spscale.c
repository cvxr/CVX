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
#include <string.h>
#include "mex.h"
#include "blksdp.h"

/* ************************************************************
   PROCEDURE realtransp - Let Y = X', where X is m x n full.
   You may also view this as going from column- towards row-stacking.
   ************************************************************ */
void realtransp(double *y, const double *x, const mwIndex m, const mwIndex n)
{
  mwIndex i, j, ji;
/* ------------------------------------------------------------
   Y(n x m), X(m x n).
   For each column y:= yj, j=0:m-1
   Let yj[i] = xji, i = 0:n-1.
   ------------------------------------------------------------ */
  for(j = 0; j < m; j++, y += n)
    for(i = 0, ji = j; i < n; i++, ji += m)
      y[i] = x[ji];            /* yij = xji */
}

/* =========================  P S D :  ========================= 
   D(d^2)x = vec(D * X * D), only at zir positions, with X sparse.
   ============================================================ */

/* ************************************************************
   PROCEDURE realDmulX - Computes Y = D*X, where X is sparse.
   INPUT
     d  - full n x n matrix D.
     xpr, xir, xjc1 - sparse vector, xjc1 nonzeros.
     first - Subscript of X(1,1), i.e. x(first:first+n^2-1) = vec(X).
     n - order of D, X.
   UPDATED
     *pxjc0 - On input, points to first nonzero with index in range
       first:first+n^2-1. On output, points to first nonzero AFTER this range.
   OUTPUT
     y  - Full n x knz output matrix. Y = D*X(:,ycols).
     ycols - length knz <= n integer array, listing all columns where X
       and hence Y is nonzero.
   RETURNS knz = length(ycols), number of nonzero columns.
   ************************************************************ */
mwIndex realdmulx(double *y, mwIndex *ycols, const double *d,
              const double *xpr, const mwIndex *xir, mwIndex *pxjc0, const mwIndex xjc1,
              const mwIndex first, const mwIndex n)
{
  mwIndex knz, jfirst, jlast, inz, i, j;
  knz = 0;        /* length(ycols) */
  jlast = 0;      /* index right after last activated column */
  y -= n;         /* point to dummy column, which will be skipped */
/* ------------------------------------------------------------
   For every nonzero xij, Let Y += D*(xij * ei*ej'), i.e.
   yj += xij * di. Let ycols list the columns where yj is nonzero.
   Store only those columns in y. (ycols is increasing)
   ------------------------------------------------------------ */
  for(inz = *pxjc0; inz < xjc1; inz++)
    if((i = xir[inz]) >= jlast){      /* we need to create new y-column */
      j = (i-first) / n;              /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all columns of X */
      ycols[knz++] = j;               /* open new column */
      y += n;
      jfirst = first + n * j;
      jlast = jfirst + n;
      i = i - jfirst;                  /* i = rem( (i-first)/n ) */
      scalarmul(y, xpr[inz], d + i * n, n);    /* yj = xij * di */
    }
    else
      addscalarmul(y, xpr[inz], d + (i-jfirst) * n, n);   /* yj += xij * di */
/* ------------------------------------------------------------
   RETURN xjc0NEW = inz, pointer to next entry in sparse x.
   knz = length(ycols), the number of nonzero columns in y.
   ------------------------------------------------------------ */
  *pxjc0 = inz;
  return knz;
}

/* ************************************************************
   PROCEDURE cpxDmulX - Computes Y = D*X, where X is sparse.
      Y,D, and X are complex matrices, stored as [vec(RE); vec(IM)].
   INPUT
     d  - full 2*(n x n) matrix D.
     xpr, xir, xjc1 - sparse vector, xjc1 nonzeros (some RE, some IM).
     first - Subscript of RE X(1,1), i.e. x(first:first+n^2-1) = vec(RE(X)),
        x(first+n^2:first+2n^2-1) = vec(IM(X)).
     n - order of D, X.
   UPDATED
     *pxjc0 - On input, points to first nonzero with index in range
       first:first+2*n^2-1. On output, first nonzero AFTER this range.
   OUTPUT
     y  - Full 2*(n x n) output matrix. y(1:n*knz) = RE D*X(:,ycols),
       y(n^2+(1:n*knz)) = IM D*X(:,ycols).
     ycols - length knz <= n integer array, listing all columns where X
       and hence Y is nonzero.
   RETURNS knz = length(ycols), number of nonzero columns.
   ************************************************************ */
mwIndex cpxdmulx(double *y, mwIndex *ycols, const double *d,
             const double *xpr, const mwIndex *xir, mwIndex *pxjc0, const mwIndex xjc1,
             const mwIndex first, const mwIndex n)
{
  mwIndex jfirst, jlast, inz, jnz, knz, i, j, icol, imgfirst, nsqr, ncols;
  double *ypr, *ypi;
  const double *dpi;
  char found;
  knz = 0;        /* length(ycols) */
  jlast = 0;      /* index right after last activated column */
  nsqr = SQR(n);
  dpi = d + nsqr;
/* ------------------------------------------------------------
   For every REAL nonzero xij, Let Y += D*(xij * ei*ej'), i.e.
   RE(yj) += xij * di, and IM(yj) += xij *IM(di).
   Let ycols list the columns where yj is nonzero due to REAL X.
   Store only those columns in y. (ycols is increasing - for now)
   ------------------------------------------------------------ */
  ypr = y-n;         /* point to dummy column, which will be skipped */
  ypi = ypr + nsqr;
  for(inz = *pxjc0; inz < xjc1; inz++)
    if((i = xir[inz]) >= jlast){      /* we need to create new y-column */
      j = (i-first) / n;              /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all real columns of X */
      ycols[knz++] = j;               /* open new column */
      ypr += n;
      ypi += n;
      jfirst = first + n * j;
      jlast = jfirst + n;
      icol = (i - jfirst) * n;               /* i = rem( (i-first)/n ) */
      scalarmul(ypr, xpr[inz], d + icol, n);               /* yj = xij * di */
      scalarmul(ypi, xpr[inz], dpi + icol, n);
    }
    else{
      icol = (i-jfirst) * n;
      addscalarmul(ypr, xpr[inz], d + icol, n); /* yj += xij * di */
      addscalarmul(ypi, xpr[inz], dpi + icol, n);
    }
/* ------------------------------------------------------------
   Finished with RE(X), yielding ncols:=knz columns in Y. Use
   jnz to browse through these existing cols while processing IMG(X).
   ------------------------------------------------------------ */
  ncols = knz;            /* #existing cols, due to RE(X) */
  jnz = 0;                /* pointer into ycols(1:ncols) */
/* ------------------------------------------------------------
   For every IMAG nonzero xij, Let Y += iD*(xij * ei*ej'), where
   iD := sqrt(-1)*D. Thus:
   RE(yj) -= xij * IM(di), and IM(yj) += xij *RE(di).
   New nonzero cols of y are listed in ycols(ncols+1:knz).
   This means that ycols consists of 2 disjoint increasing lists.
   ------------------------------------------------------------ */
  imgfirst = first + nsqr;
  for(; inz < xjc1; inz++)
    if((i = xir[inz]) >= jlast){      /* we need to create new y-column */
      j = (i-imgfirst) / n;           /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all columns of X */
      jfirst = imgfirst + n * j;
      jlast = jfirst + n;
      icol = (i - jfirst) * n;               /* i = rem( (i-first)/n ) */
/* ------------------------------------------------------------
   1st time to use column j while processing IM(X). See whether
   yj already in ycols (=>found=1). Otherwise, create new yj.
   ------------------------------------------------------------ */
      for(found = 0; jnz < ncols; jnz++)
        if(ycols[jnz] >= j){
          found = (ycols[jnz] == j);
          break;                        /* ycols(1:ncols) is increasing */
        }
      if(!found){
        ypr = y + n * knz;              /* IM(X) results in new yj */
        ypi = ypr + nsqr;
        ycols[knz++] = j;               /* open new column */
        scalarmul(ypr, -xpr[inz], dpi + icol, n);     /* yj = xij * (i*di) */
        scalarmul(ypi, xpr[inz], d + icol, n);
      }
      else{
        ypr = y + n * jnz;             /* yj already exists */
        ypi = ypr + nsqr;
        addscalarmul(ypr, -xpr[inz], dpi + icol, n); /* yj += xij * (i*di) */
        addscalarmul(ypi, xpr[inz], d + icol, n);
      }
    }
    else{
      icol = (i-jfirst) * n;
      addscalarmul(ypr, -xpr[inz], dpi + icol, n); /* yj += xij * (i*di) */
      addscalarmul(ypi, xpr[inz], d + icol, n);
    }
/* ------------------------------------------------------------
   RETURN xjc0NEW = inz, pointer to next entry in sparse x.
   knz = length(ycols), the number of nonzero columns in y.
   ------------------------------------------------------------ */
  *pxjc0 = inz;
  return knz;
}

/* ************************************************************
   PROCEDURE sprealdxd - Computes Z = D*sym(X)*D with D = D' and
        sym(X) = (X+X')/2. X sparse, and Z has pre-chosen sparsity
        structure zir (typically subset of tril-entries).
   INPUT
     zir, zjc0, zjc1 - Target sparsity structure. We compute z(i) only for
       i in zir(zjc0:zjc1-1) with i < first + n^2.
       zjc0 such that zir(zjc0) >= first.
     d      - full n x n matrix D.
     xpr, xir, xjc1 - sparse vector, xjc1 nonzeros.
     first  - Subscript of X(1,1) and Z(1,1), i.e.
       x(first:first+n^2-1) = vec(X).
     n      - order of D, X.
   UPDATED
     *pxjc0 - On input, points to first nonzero with index in range
       first:first+n^2-1. On output, points to first nonzero AFTER this range.
   OUTPUT
     Z - lenfull vector. Only z(first:first+n^2-1) is affected, and from this
       only the indices listed in zir. (Remaining entries unaffected.)
   WORK
     dxcols - length n integer array.
     fwork  - 2 * n^2 vector.
   ************************************************************ */
void sprealdxd(double *z, const mwIndex *zir, const mwIndex znnz,
               const double *d,
               const double *xpr, const mwIndex *xir, mwIndex *pxjc0, const mwIndex xjc1,
               const mwIndex first, const mwIndex n, double *fwork, mwIndex *dxcols)
{
  mwIndex inz, i, icol, j, jfirst, jlast, m;
  double *xd;
  const double *dj, *xdj;
/* ------------------------------------------------------------
   Partition 2*n^2 WORKING array fwork into [fwork(n^2), xd(n^2)].
   ------------------------------------------------------------ */
  xd = fwork + SQR(n);
/* ------------------------------------------------------------
   Let dxcols={j: X(:,j)!= all-0}, m = |dxcols|,
   fwork = D*X(:,dxcols).
   ------------------------------------------------------------ */
  m = realdmulx(fwork, dxcols, d, xpr, xir, pxjc0, xjc1, first, n);
/* ------------------------------------------------------------
   Let xd = fwork' = X(:,dxcols)'*D.
   ------------------------------------------------------------ */
  realtransp(xd, fwork, n,m);
/* ------------------------------------------------------------
   Let fwork = D(dxcols,:).     Both xd and fwork are m x n.
   ------------------------------------------------------------ */
  for(j = 0, inz = 0; j < n; j++, d += n)
    for(i = 0; i < m; i++)
      fwork[inz++] = d[dxcols[i]];        /* fwork(i,j) = d(dxcols(i),j) */
/* ------------------------------------------------------------
   For all target subscripts (i,j), let
   zij = (D*sym(X)*D)_ij = [ DXD_ij + DXD_ji ] /2.
   Note that DXD_ij = xd(:,i)' * fwork(:,j).   (m mults)
   ------------------------------------------------------------ */
  jlast = 0;      /* index right after last activated column */
  for(inz = 0; inz < znnz; inz++){
    if((i = zir[inz]) >= jlast){      /* move to new z-column */
      j = (i-first) / n;              /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all columns of Z */
      jfirst = first + n * j;
      jlast = jfirst + n;
      icol = m * j;
      dj = fwork + icol;
      xdj = xd + icol;
      if(i - jfirst == j)             /* zjj = xdj' * dj */
        z[i] = realdot(xdj,dj, m);
      else{              /* zij = (xdj'*di + xdi'*dj) / 2 */
        icol = (i - jfirst) * m;
        z[i] = (realdot(xdj,fwork + icol, m) + realdot(xd + icol, dj, m)) / 2;
      }
    }
    else{               /* zij = (xdj'*di + xdi'*dj) / 2 */
      icol = (i - jfirst) * m;
      z[i] = (realdot(xdj,fwork + icol, m) + realdot(xd + icol, dj, m)) / 2;
    }
  }
}

/* ************************************************************
   PROCEDURE spcpxdxd - Computes Z = D*herm(X)*D with D = D' and
        herm(X) = (X+X')/2. X sparse, and Z has pre-chosen sparsity
        structure zir (typically subset of tril-entries).
      Z,D, and X are complex matrices, stored as [vec(RE); vec(IM)].
   INPUT
     zir, zjc0, zjc1 - Target sparsity structure. We compute z(i) only for
       i in zir(zjc0:zjc1-1) with i < first + 2*n^2.
       zjc0 such that zir(zjc0) >= first.
     d      - full 2*(n x n) matrix D.
     xpr, xir, xjc1 - sparse vector, xjc1 nonzeros.
     first  - Subscript of X(1,1) and Z(1,1), i.e.
       x(first:first+n^2-1) = vec(real(X)), x(first+n^2:end)=vec(imag(X)).
     n      - order of D, X.
   UPDATED
     *pxjc0 - On input, points to first nonzero with index in range
       first:first+2*n^2-1. On output, points just BEYOND.
   OUTPUT
     Z - lenfull vector. Only z(first:first+2*n^2-1) is affected, and from this
       only the indices listed in zir. (Other entries unaffected.)
   WORK
     dxcols - length n integer array.
     fwork  - 4 * n^2 vector.
   RETURNS zjc0_NEW
   ************************************************************ */
void spcpxdxd(double *z, const mwIndex *zir, const mwIndex znnz,
              const double *d,
              const double *xpr, const mwIndex *xir, mwIndex *pxjc0, const mwIndex xjc1,
              const mwIndex first, const mwIndex n, double *fwork, mwIndex *dxcols)
{
  mwIndex inz, i, icol, j, jfirst, jlast, m, nsqr, imgfirst;
  double *dx, *dxpi, *fworkpi;
  double zi;
  const double *dj, *djpi, *djx, *djxpi;
  nsqr = SQR(n);
/* ------------------------------------------------------------
   Partition 4*n^2 WORKING array fwork into [fwork(2*n^2), dxRows(2*n^2)].
   ------------------------------------------------------------ */
  fworkpi = fwork + nsqr;
  dx = fworkpi + nsqr;
  dxpi = dx + nsqr;
/* ------------------------------------------------------------
   Let dxcols={j: X(:,j)!= all-0}, m = |dxcols|,
   fwork = RE(D*X(:,dxcols)), fwork+n^2 = IM D*X(:,dxcols).
   ------------------------------------------------------------ */
  m = cpxdmulx(fwork, dxcols, d, xpr, xir, pxjc0, xjc1, first, n);
/* ------------------------------------------------------------
   Let dxRows = D*X(:,dxcols) in row-stacking (instead of usual col-stack).
   ------------------------------------------------------------ */
  realtransp(dx, fwork, n,m);
  realtransp(dxpi, fworkpi, n,m);
/* ------------------------------------------------------------
   Let fwork = D(dxcols,:) in (usual) column stacking; m x n.
   ------------------------------------------------------------ */
  for(j = 0, inz = 0; j < n; j++, d += n)
    for(i = 0; i < m; i++)
      fwork[inz++] = d[dxcols[i]];        /* fwork(i,j) = RE d(dxcols(i),j) */
  for(j = 0, inz = 0; j < n; j++, d += n)
    for(i = 0; i < m; i++)
      fworkpi[inz++] = d[dxcols[i]];     /* fworkpi(i,j) = IM d(dxcols(i),j) */
/* ------------------------------------------------------------
   For all target subscripts (i,j) in RE Z, let
   zij = RE (D*herm(X)*D)_ij = RE [ DXD_ij + DXD_ji ] /2.
   Note that DXD_ij = dx(i,:) * fwork(:,j).   (m mults)
   Since dx is in row-stacking, we can simply browse thru memory.
   ------------------------------------------------------------ */
  jlast = 0;      /* index right after last activated column */
  for(inz = 0; inz < znnz; inz++){
    if((i = zir[inz]) >= jlast){      /* move to new z-column */
      j = (i-first) / n;              /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all columns of RE(Z) */
      jfirst = first + n * j;
      jlast = jfirst + n;
      icol = m * j;
      dj = fwork + icol;
      djpi = dj + nsqr;
      djx = dx + icol;
      djxpi = djx + nsqr;
      if(i - jfirst == j)             /* zjj = djx * dj - IMdjx*IMdj * */
        z[i] = realdot(djx,dj, m) - realdot(djxpi,djpi, m);
      else{              /* zij = RE (dix*dj + conj(djx*di)) / 2 */
        icol = (i - jfirst) * m;
        zi = realdot(djx,fwork + icol, m) - realdot(djxpi,fworkpi+icol,m);
        zi += realdot(dx + icol, dj, m) - realdot(dxpi+icol, djpi, m);
        z[i] = zi / 2;
      }
    }
    else{               /* zij = RE (dix*dj + conj(djx*di)) / 2 */
      icol = (i - jfirst) * m;
      zi = realdot(djx,fwork + icol, m) - realdot(djxpi,fworkpi+icol,m);
      zi += realdot(dx + icol, dj, m) - realdot(dxpi+icol, djpi, m);
      z[i] = zi / 2;
    }
  }
/* ------------------------------------------------------------
   Remains: target subscripts (i,j) in IM Z.
   zij = IM (D*herm(X)*D)_ij = IM [ DXD_ij + DXD_ji ] /2.
   Note that DXD_ij = dx(i,:) * fwork(:,j).   (m mults)
   Since dx is in row-stacking, we can simply browse thru memory.
   ------------------------------------------------------------ */
  imgfirst = first + nsqr;
  for(; inz < znnz; inz++){
    if((i = zir[inz]) >= jlast){      /* move to new z-column */
      j = (i-imgfirst) / n;              /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all columns of IM(Z) */
      jfirst = imgfirst + n * j;
      jlast = jfirst + n;
      icol = m * j;
      dj = fwork + icol;
      djpi = dj + nsqr;
      djx = dx + icol;
      djxpi = djx + nsqr;
      mxAssert(i - jfirst != j,"");          /* Herm => diag(imag(Z))=0 */
      icol = (i - jfirst) * m;     /* zij = IM (dix*dj + conj(djx*di)) / 2 */
      zi = realdot(dx + icol, djpi, m) + realdot(dxpi+icol, dj, m);
      zi -= realdot(djx,fworkpi + icol, m) + realdot(djxpi,fwork+icol,m);
      z[i] = zi / 2;
    }
    else{               /* zij = IM (dix*dj + conj(djx*di)) / 2 */
      icol = (i - jfirst) * m;
      zi = realdot(dx + icol, djpi, m) + realdot(dxpi+icol, dj, m);
      zi -= realdot(djx,fworkpi + icol, m) + realdot(djxpi,fwork+icol,m);
      z[i] = zi / 2;
    }
  }
}

/* ************************************************************
   PROCEDURE spsqrscale - Computes z = D(d^2)x for PSD-part.
   INPUT
     zir, zjc1 - Target sparsity structure. We compute z(i) only for
       i in zir(0:zjc1-1).
     d      - lenud vector containing full n x n matrices Dk (=Ud'*Ud).
     xpr, xir, xjc1 - sparse vector, xjc1 nonzeros.
     xjc0 - Points to first nonzero with index in range
       blkstart[0]:blkstart[psdN].
     blkstart  - length psdN+1 array such that 
       x(blkstart[k]:blkstart[k+1]) = vec(Xk) for all PSD blocks k.
     xblk - length psdDim array, with k = xblk(i-blkstart[0]) iff
       blkstart[k] <= i < blkstart[k+1], k=0:nblk-1.
       psdDim:=blkstart[end]-blkstart[0].
     psdNL     - length psdN array, is K.s.
     rpsdN     - number of real symmetric PSD blocks.
   OUTPUT
     z - lenfull vector. Only indices listed in zir are affected.
       Other entries are unaffected (not even set to 0).
     blks - length blksnnz <= length(K.s), lists PSD-blocks where
       z has nonzeros.
   WORK
     fwork - 2 * (rMaxn^2 + 2*hMaxn^2) vector.  (<= 4*max(K.s)^2)
     iwork - length n=MAX(rMaxn,hMaxn) integer array. (n=max(K.s))
   RETURNS blksnnz, number of entries written into blks.
   ************************************************************ */
mwIndex spsqrscale(double *z, mwIndex *blks, const mwIndex *zjc, const mwIndex *zir,
               const mwIndex *znnz, const double *d,
               const mwIndex *xir, const double *xpr, mwIndex xjc0, const mwIndex xjc1,
               const mwIndex *blkstart, const mwIndex *xblk, const mwIndex *psdNL,
               const mwIndex rpsdN, double *fwork, mwIndex *iwork)
{
  mwIndex blksnnz, k;

  blksnnz = 0;
  d -= blkstart[0];        /* Now, d + blkstart[k] gives Dk */
  while(xjc0 < xjc1){
/* ------------------------------------------------------------
   For each nonzero block xblk[.] in x, Let zk = vec(tril(Dk * Xk * Dk)).
   ------------------------------------------------------------ */
    if((k = xblk[xir[xjc0]]) >= rpsdN)
      break;                                /* 1st nz Hermitian PSD block */
    sprealdxd(z, zir+zjc[k], znnz[k], d + blkstart[k], xpr, xir, &xjc0, xjc1,
              blkstart[k], psdNL[k], fwork, iwork);
    blks[blksnnz++] = k;
  }
/* ------------------------------------------------------------
   Hermitian PSD blocks
   ------------------------------------------------------------ */
  while(xjc0 < xjc1){
    k = xblk[xir[xjc0]];
    spcpxdxd(z, zir+zjc[k], znnz[k], d + blkstart[k], xpr, xir, &xjc0, xjc1,
             blkstart[k], psdNL[k], fwork, iwork);
    blks[blksnnz++] = k;
  }
  return blksnnz;
}

#ifdef OLDSEDUMI
/* =========================  P S D :  ========================= 
   D(d)x = vec(U' * X * U), only at zir positions.
   Typically, zir is subset of tril (or of triu).
   ============================================================ */

/* ************************************************************
   PROCEDURE sprealutxu - Computes Z(perm,perm) = U' * X * U.
        Z has pre-chosen sparsity structure zir
	(typically subset of tril-entries [or of triu-entries]).
   INPUT
     zir, zjc0, zjc1 - Target sparsity structure. We compute z(i) only for
       i in zir(zjc0:zjc1-1) with i < first + n^2.
       zjc0 such that zir(zjc0) >= first.
     u      - full n x n matrix, triu(u) = U.
     invperm - length n array, invperm(perm)=0:n-1. Used to translate Z-index
       to u- and x-index.
     x - n x n full symmetric matrix containing triu(X).
     first  - Subscript of X(1,1) and Z(1,1), i.e.
       x(first:first+n^2-1) = vec(X).
     n      - order of U, X.
   OUTPUT
     Z - lenfull vector. Only z(first:first+n^2-1) is affected, and from this
       only the indices listed in zir. (Other entries unaffected.)
   WORK
     fwork  - 2 * n^2 vector.
   RETURNS zjc0_NEW
   ************************************************************ */
mwIndex sprealutxu(double *z, const mwIndex *zir, const mwIndex zjc0, const mwIndex zjc1,
               const double *u, const mwIndex *invperm, const double *x,
               const mwIndex first, const mwIndex n, double *fwork)
{
  mwIndex inz, i, irow, j, icol, jcol, jfirst, jlast, m;
  double *xu;
  const double *uj, *xuj;
/* ------------------------------------------------------------
   Partition 2*n^2 WORKING array fwork into [fwork(n^2), xu(n^2)].
   ------------------------------------------------------------ */
  xu = fwork + SQR(n);
/* ------------------------------------------------------------
   Let fwork = X in full symmetric storage
   ------------------------------------------------------------ */
  memcpy(fwork, x, SQR(n) * sizeof(double));
  triu2sym(fwork,n);
/* ------------------------------------------------------------
   Compute xu = triu(fwork'*U):   1^2 + 2^2 + ... + (n-1)^2 + n^2 mults.
   For each j, let uj point to u(0,j); m is
   remaining length of column j, i.e. m=j.      (SEE realutxu.)
   ------------------------------------------------------------ */
  inz = 0; jcol = 0;
  for(j = 0, m = n; j < n; j++, jcol += n){
    for(i = 0, icol = 0; i <= j; i++, icol += n)
      xu[inz++] = realdot(fwork+icol,u+jcol,j+1);  /* x(:,i)'*u(0:j,j) */
    inz += --m;                                    /* skip tril */
  }
/* ------------------------------------------------------------
   For all target subscripts (i,j), let
   zij = (U'*XU)_ij [ = (U'*XU)_ji ]
   ------------------------------------------------------------ */
  jlast = 0;                  /* index right after last activated column */
  for(inz = zjc0; inz < zjc1; inz++){
    if((i = zir[inz]) >= jlast){      /* move to new z-column */
      j = (i-first) / n;              /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all columns of Z */
      jfirst = first + n * j;
      jlast = jfirst + n;
      j = invperm[j];                 /* change to U- and X-index */
      jcol = n * j;
      uj = u + jcol;
      xuj = xu + jcol;
    }
    irow = invperm[i - jfirst];       /* U- and X-index */
    if(irow <= j)            /* zij = u(0:i,i)' * xu(0:i,j) IF i<=j */
      z[i] = realdot(u + irow * n, xuj, irow + 1);
    else{                    /* zij = u(0:j,j)' * xu(0:j,i) IF i>j */
      z[i] = realdot(uj, xu + irow * n, j + 1);
    }
  }
/* ------------------------------------------------------------
   RETURN inz, next position in zir (after range first:first+n^2-1).
   ------------------------------------------------------------ */
  return inz;
}

/* ************************************************************
   PROCEDURE spcpxutxu - Computes Z(perm,perm) = U' * X * U.
        Z has pre-chosen sparsity structure zir
	(typically subset of tril-entries [or of triu-entries]).
        U,X and Z are complex matrices, stored as [vec(RE); vec(IM)].
   INPUT
     zir, zjc0, zjc1 - Target sparsity structure. We compute z(i) only for
       i in zir(zjc0:zjc1-1) with i < first + 2*n^2.
       zjc0 such that zir(zjc0) >= first.
     u      - full 2*(n x n) matrix, triu(u) = U.
     invperm - length n array, invperm(perm)=0:n-1. Used to translate Z-index
       to u- and x-index.
     x - 2*(n x n) full symmetric matrix containing triu(X).
     first  - Subscript of X(1,1) and Z(1,1), i.e.
       x(first:first+n^2-1) = vec(RE X), x(first+n^2:end) = vec(IM(X))..
     n      - order of U, X.
   OUTPUT
     Z - lenfull vector. Only z(first:first+2n^2-1) is affected, and from this
       only the indices listed in zir. (Other entries unaffected.)
   WORK
     fwork  - 4 * n^2 vector.
   RETURNS zjc0_NEW
   ************************************************************ */
mwIndex spcpxutxu(double *z, const mwIndex *zir, const mwIndex zjc0, const mwIndex zjc1,
              const double *u, const mwIndex *invperm, const double *x,
              const mwIndex first, const mwIndex n, double *fwork)
{
  mwIndex inz, i, irow, j, icol, jcol, jfirst, jlast, m, nsqr, imgfirst;
  double *xu, *xupi, *fworkpi;
  const double *uj, *ujpi, *xuj, *xujpi, *upi;
/* ------------------------------------------------------------
   Partition 4*n^2 WORKING array fwork into [fwork(2*n^2), xu(2*n^2)].
   ------------------------------------------------------------ */
  nsqr = SQR(n);
  upi = u + nsqr;
  fworkpi = fwork + nsqr;
  xu = fworkpi + nsqr;
  xupi = xu + nsqr;
/* ------------------------------------------------------------
   Let fwork = X in full hermitian storage
   ------------------------------------------------------------ */
  memcpy(fwork, x, 2*nsqr * sizeof(double));
  triu2herm(fwork,fworkpi,n);
/* ------------------------------------------------------------
   Compute xu = triu(fwork'*U)   ( =triu(XU), since X is Hermitian.)   
   For each j, let uj point to u(0,j); m is
   remaining length of column j, i.e. m=j.      (SEE prpiutxu.)
   ------------------------------------------------------------ */
  inz = 0; jcol = 0;
  for(j = 0, m = n; j < n; j++, jcol += n){
    for(i = 0, icol = 0; i <= j; i++, icol += n){
/* xu(i,j) = x(:,i)'*u(0:j,j), assume IM u(j,j) == 0 */
      xu[inz] = realdot(fwork+icol,u+jcol,j+1) +
        realdot(fworkpi+icol,upi+jcol,j);
      xupi[inz] = realdot(fwork+icol,upi+jcol,j) -
        realdot(fworkpi+icol,u+jcol,j+1);
      inz++;
    }
    inz += --m;                                /* skip tril */
  }
/* ------------------------------------------------------------
   For all target subscripts (i,j) in RE Z, let
   zij = RE (U'*XU)_ij [ = RE (U'*XU)_ji ]
   ------------------------------------------------------------ */
  jlast = 0;                  /* index right after last activated column */
  for(inz = zjc0; inz < zjc1; inz++){
    if((i = zir[inz]) >= jlast){      /* move to new z-column */
      j = (i-first) / n;              /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all columns of RE Z */
      jfirst = first + n * j;
      jlast = jfirst + n;
      j = invperm[j];                 /* change to U- and X-index */
      jcol = n * j;
      uj = u + jcol;
      ujpi = uj + nsqr;
      xuj = xu + jcol;
      xujpi = xuj + nsqr;
    }
    irow = invperm[i - jfirst];       /* U- and X-index */
    icol = irow * n;
    if(irow <= j)            /* zij = RE u(0:i,i)' * xu(0:i,j) IF i<=j */
      z[i] = realdot(u+icol, xuj, irow+1) + realdot(upi+icol, xujpi, irow);
    else{                    /* zij = RE u(0:j,j)' * xu(0:j,i) IF i>j */
      z[i] = realdot(uj, xu+icol, j+1) + realdot(ujpi, xupi+icol, j);
    }
  }
/* ------------------------------------------------------------
   Remains target subscripts (i,j) in IM Z. Let
   zij = IM (U'*XU)_ij [ = -IM (U'*XU)_ji ]
   ------------------------------------------------------------ */
  imgfirst = first + nsqr;
  for(; inz < zjc1; inz++){
    if((i = zir[inz]) >= jlast){      /* move to new z-column */
      j = (i-imgfirst) / n;              /* col j = floor( (i-first)/n ) */
      if(j >= n)
        break;                        /* finished all columns of IM Z */
      jfirst = imgfirst + n * j;
      jlast = jfirst + n;
      j = invperm[j];                 /* change to U- and X-index */
      jcol = n * j;
      uj = u + jcol;
      ujpi = uj + nsqr;
      xuj = xu + jcol;
      xujpi = xuj + nsqr;
    }
    irow = invperm[i - jfirst];       /* U- and X-index */
    icol = irow * n;
    if(irow <= j)            /* zij = IM u(0:i,i)' * xu(0:i,j) IF i<=j */
      z[i] = realdot(u+icol, xujpi, irow+1) - realdot(upi+icol, xuj, irow);
    else                     /* zij = -IM u(0:j,j)' * xu(0:j,i) IF i>j */
      z[i] = realdot(ujpi, xu+icol, j) - realdot(uj, xupi+icol, j+1);
  }
/* ------------------------------------------------------------
   RETURN inz, next position in zir (after range first:first+2*n^2-1).
   ------------------------------------------------------------ */
  return inz;
}

/* ************************************************************
   PROCEDURE blkpsdscale - Computes z = D(d)x = U'*X*U for PSD-part.
     Typically, X = D(d)aj = vec(triu(U*Aj*U')).
   INPUT
     zir, zjc1 - Target sparsity structure. We compute z(i) only for
       i in zir(0:zjc1-1).
     u      - lenud vector containing full n x n matrices Uk
     invperm - lenpsd array, containing inverse perm for Uk for each psd-block
     x      - vector of full blocks, viz. those listed in xblk.
     xblk, blkjc0, blkjc1 - the psd-blocks stored in x.
     blkstart  - length psdN+1 array such that
       z(blkstart[k]:blkstart[k+1]) = vec(Zk) for all PSD blocks k.
     psdNL     - length psdN array, is K.s.
     cumpsdNL  - cumsum([0, psdNL]), i.e. start of block k in invperm.
     rpsdN     - number of real symmetric PSD blocks.
   OUTPUT (PRE-INITIALIZED)
     z - lenfull vector. Only indices listed in zir are affected.
       Other entries are unaffected (not even set to 0).
   WORK
     fwork  - 2 * (rMaxn^2 + 2*hMaxn^2) vector.  (<= 4*max(K.s)^2)
   RETURNS zjc0_NEW (should be zjc1)
   ************************************************************ */
mwIndex blkpsdscale(double *z, const mwIndex *zir, const mwIndex zjc1,
		const double *u, const mwIndex *invperm, const double *x,
		const mwIndex *xblk, const mwIndex blkjc0, const mwIndex blkjc1,
		const mwIndex *blkstart, const mwIndex *psdNL, const mwIndex *cumpsdNL,
		const mwIndex rpsdN, double *fwork)
{
  mwIndex inz, knz, k, nk, lastzsub;

  inz = 0;                 /* points into zir */
  if(inz >= zjc1)
    return inz;            /* return if no target subscripts */
  lastzsub = zir[zjc1-1];
  u -= blkstart[0];        /* Now, u + blkstart[k] gives Uk */
  for(knz = blkjc0; knz < blkjc1; knz++){
/* ------------------------------------------------------------
   For each nonzero block xblk[k] in x, Let zk = vec(tril(Dk * Xk * Dk)).
   ------------------------------------------------------------ */
    k = xblk[knz];
    if(k >= rpsdN)
      break;                                /* 1st nz Hermitian PSD block */
    if(lastzsub < blkstart[k])
      break;                                /* already beyond last target i*/
    nk = psdNL[k];
    while(zir[inz] < blkstart[k])           /* move to block k, viz. */
      z[zir[inz++]] = 0.0;                  /* all-0 where x is all-0 */
    inz = sprealutxu(z, zir, inz, zjc1, u + blkstart[k], invperm+cumpsdNL[k],
		     x, blkstart[k], nk, fwork);
    x += SQR(nk);                           /* point to next nonzero block */
  }
/* ------------------------------------------------------------
   Hermitian PSD blocks
   ------------------------------------------------------------ */
  for(; knz < blkjc1; knz++){
    k = xblk[knz];
    if(lastzsub < blkstart[k])
      break;                                /* already beyond last target i*/
    nk = psdNL[k];
    while(zir[inz] < blkstart[k])           /* move to block k, viz. */
      z[zir[inz++]] = 0.0;                  /* all-0 where x is all-0 */
    inz = spcpxutxu(z, zir, inz, zjc1, u + blkstart[k], invperm+cumpsdNL[k],
                    x, blkstart[k], nk, fwork);
    x += 2*SQR(nk);                         /* point to next nonzero block */
  }
/* ------------------------------------------------------------
   No remaining x blocks ==> remaining z-entries all-0.
   ------------------------------------------------------------ */
  for(; inz < zjc1; inz++)
    z[zir[inz]] = 0.0;                  /* all-0 where x is all-0 */
/* ------------------------------------------------------------
   RETURN next zir-position (should be zjc1)
   ------------------------------------------------------------ */
  return inz;       /* arrived at end of z-vector */
}
#endif
