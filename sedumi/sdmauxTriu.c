/* ************************************************************
   MODULE sdmaux*.c  -- Several low-level subroutines for the
   mex-files in the Self-Dual-Minimization package.

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
#include "blksdp.h"

/* ============================================================
   TRIU/TRIL RELATED
   ============================================================ */
/* ************************************************************
   PROCEDURE tril2sym -- assign R(i,j) = R(j,i) for all j>i.
   INPUT n - order of n x n matrix R.
   UPDATED r - on output, the strict lower triangular is copied
       to the strict upper triangular.
   ************************************************************ */
void tril2sym(double *r, const mwIndex n)
{
 mwIndex colp,i,j;

 /* ------------------------------------------------------------
    r points to R(:,i);     r+colp = R(:,j).
    ------------------------------------------------------------ */
 for(i=0; i<n; r += n, i++)
   for(colp = n + i, j=i+1; j<n; j++, colp += n)
     r[colp] = r[j];                          /* R(i,j) = R(j,i) */
}

/* ************************************************************
   PROCEDURE tril2herm -- Given R = [RE R, IM R],
     assign RE R(i,j) = RE R(j,i) and IM R(i,j) = - IM R(j,i) for all j>i.
   INPUT n - order of 2*(n x n) matrix R.
   UPDATED r,rpi - on output, is made Hermitian (sym, skew-sym resp).
   ************************************************************ */
void tril2herm(double *r, double *rpi, const mwIndex n)
{
  mwIndex colp,i,j;
/* ------------------------------------------------------------
   First, make the real block symmetric. Then, make the imaginary
   part skew-symmetric.
    ------------------------------------------------------------ */
  tril2sym(r,n);
/* ------------------------------------------------------------
   r points to R(:,i);     r+colp = R(:,j).
   ------------------------------------------------------------ */
  for(i = 0; i < n; rpi += n, i++){
    rpi[i] = 0.0;                                 /* zero-diagonal */
    for(colp = n + i, j = i + 1; j < n; j++, colp += n)
      rpi[colp] = -rpi[j];                        /* R(i,j) = -R(j,i) */
  }
}

/* ************************************************************
   PROCEDURE triu2sym -- assign R(j,i) = R(i,j) for all j>i.
   INPUT n - order of n x n matrix R.
   UPDATED r - on output, the strict lower triangular is copied
       to the strict upper triangular.
   ************************************************************ */
void triu2sym(double *r, const mwIndex n)
{
 mwIndex colp,i,j;

 /* ------------------------------------------------------------
    r points to R(:,i);     r+colp = R(:,j).
    ------------------------------------------------------------ */
 for(i = 0; i < n; r += n, i++)
   for(colp = n + i, j=i+1; j<n; j++, colp += n)
     r[j] = r[colp];                          /* R(j,i) = R(i,j) */
}

/* ************************************************************
   PROCEDURE triu2herm -- Given R = [RE R, IM R],
     assign RE R(i,j) = RE R(j,i) and IM R(i,j) = - IM R(j,i) for all j<i.
   INPUT n - order of 2*(n x n) matrix R.
   UPDATED r, rpi - on output, is made Hermitian (sym and skewsym resp).
   ************************************************************ */
void triu2herm(double *r, double *rpi, const mwIndex n)
{
 mwIndex colp,i,j;
 /* ------------------------------------------------------------
    First, make the real block symmetric. Then, let r point to
    the imaginary part and make that skew-symmetric.
    ------------------------------------------------------------ */
 triu2sym(r,n);
/* ------------------------------------------------------------
   rpi points to R(:,i);     rpi+colp = R(:,j).
   ------------------------------------------------------------ */
 for(i = 0; i < n; rpi += n, i++){
   rpi[i] = 0.0;                                 /* zero-diagonal */
   for(colp = n + i, j=i+1; j < n; j++, colp += n)
     rpi[j] = -rpi[colp];                          /* R(j,i) = -R(i,j) */
 }
}

/* ************************************************************
   PROCEDURE uperm - Let triu(Y) = triu(U(:,perm)),
                     leaving tril(Y,-1) undefined.
   INPUT
     u - nxn input matrix
     perm - length n pivot ordering
     n -   order
   OUTPUT
     y - triu(y) = triu(u(:,perm).
   ************************************************************ */
void uperm(double *y, const double *u, const mwIndex *perm, const mwIndex n)
{
  mwIndex j;
  const double *uj;

  for(j = 0; j < n; y += n){
    uj = u + perm[j] * n;
    memcpy(y, uj, (++j) * sizeof(double));
  }
}
