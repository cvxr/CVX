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
#include "triuaux.h"
#include "blksdp.h"

/* ************************************************************
   PROCEDURE matperm - Let Y = X(perm,perm)
   INPUT
     x -  n x n input matrix
     perm - length n permutation
     n - order
   OUTPUT
     y - n x n matrix, Y = X(perm,perm)
   ************************************************************ */
void matperm(double *y, const double *x, const mwIndex *perm, const mwIndex n)
{
  mwIndex i,j,inz;
  const double *xj;

  inz = 0;
  for(j = 0; j < n; j++){
    xj = x + perm[j] * n;
    for(i = 0; i < n; i++)
      y[inz++] = xj[perm[i]];
  }
}

/* ************************************************************
   PROCEDURE invmatperm - Let Y(perm,perm) = X
   ************************************************************ */
void invmatperm(double *y, const double *x, const mwIndex *perm, const mwIndex n)
{
  mwIndex i,j,inz;
  double *yj;

  inz = 0;
  for(j = 0; j < n; j++){
    yj = y + perm[j] * n;
    for(i = 0; i < n; i++)
      yj[perm[i]] = x[inz++];
  }
}

/* ************************************************************
   PROCEDURE realltxl: Computes L'*X*L in only n*(n+1)^2/2
   multiplications.
   INPUT
     l - n x n full matrix, tril(l) = L.
     x - n x n full symmetric matrix: instead of X*L, we actually
       compute X'*L.
     n - order.
   OUTPUT
     y - n x n full matrix; tril(Y) = tril(L'*X*L).
   WORK
     xl - length n^2 working vector, to store tril(X'*L).
   ************************************************************ */
void realltxl(double *y, const double *ld, const double *x,
              const mwIndex n, double *xld)
{
  mwIndex inz, i,j, m, icol, jcol;
/* ------------------------------------------------------------
   Compute xl = tril(X'*L):   n^2 + (n-1)^2 + .. + 1^2 mults.
   For each j, let xj and lj point to x(j,j) and l(j,j); m is
   remaining length of column j, i.e. m=n-j.
   ------------------------------------------------------------ */
  inz = 0; jcol = 0;
  for(j = 0, m = n; j < n; j++, m--, jcol += n+1){
    inz += j;                                 /* go to tril */
    for(i = j, icol = jcol; i < n; i++, icol += n)
      xld[inz++] = realdot(x+icol,ld+jcol,m);   /*xl(i,j) = x(:,i)'*l(j:n,j) */
  }
/* ------------------------------------------------------------
   Compute y = tril(L'*XL):  1*n + 2*(n-1) + ... + n*1 mults.
   ------------------------------------------------------------ */
  inz = 0; jcol = 0;
  for(j = 0; j < n; j++, jcol += n+1){
    inz += j;                                 /* go to tril */
    for(i = j, icol = jcol; i < n; i++, icol += n+1){
      y[inz] = realdot(ld + icol,xld+inz,n-i); /* yij = l(i:n,i)' * xl(:,j) */
      inz++;
    }
  }
}

/* ************************************************************
   PROCEDURE prpiltxl: Computes L'*X*L, with X,L complex, diag(L) real.
   INPUT
     l - n x n full matrix, tril(l) = L. Assume IM diag(L) is all-0.
     x - n x n full Hermitian matrix: instead of X*L, we actually
       compute X'*L.
     n - order.
   OUTPUT
     y - n x n full matrix; tril(Y) = tril(L'*X*L).
   WORK
     xld - length 2 * n^2 working vector, to store tril(X'*L).
   ************************************************************ */
void prpiltxl(double *y,double *ypi, const double *ld,const double *ldpi,
              const double *x,const double *xpi,
              const mwIndex n, double *xld)
{
  mwIndex inz, i,j, m, icol, jcol;
  double *xldpi;
/* ------------------------------------------------------------
   Partition working array xld into real and imaginary part
   ------------------------------------------------------------ */
  xldpi = xld + SQR(n);
/* ------------------------------------------------------------
   Compute xl = tril(X'*L)        (= tril(XL), since X is Hermitian.)
   For each j, let xj and lj point to x(j,j) and l(j,j); m is
   remaining length of column j, i.e. m=n-j.
   ------------------------------------------------------------ */
  inz = 0; jcol = 0;
  for(j = 0, m = n; j < n; j++, m--, jcol += n+1){
    inz += j;                                 /* go to tril */
    for(i = j, icol = jcol; i < n; i++, icol += n){
/* xl(i,j) = x(:,i)'*l(j:n,j), assume IM l(j,j) == 0 */
      xld[inz] = realdot(x+icol,ld+jcol,m)
        + realdot(xpi+icol+1,ldpi+jcol+1,m-1);
      xldpi[inz] = realdot(x+icol+1,ldpi+jcol+1,m-1)
        - realdot(xpi+icol,ld+jcol,m);
      inz++;
    }
  }
/* ------------------------------------------------------------
   Compute y = tril(L'*XL)
   ------------------------------------------------------------ */
  inz = 0; jcol = 0;
  for(j = 0; j < n; j++, jcol += n+1){
    inz += j;                                 /* go to tril */
    for(i = j, icol = jcol; i < n; i++, icol += n+1){
/* yij = l(i:n,i)' * xl(:,j). Assume IM l(i,i) == 0 */
      y[inz] = realdot(ld + icol,xld+inz,n-i)
        + realdot(ldpi + icol+1,xldpi+inz+1,n-i-1);
      ypi[inz] = realdot(ld + icol,xldpi+inz,n-i)
        - realdot(ldpi+icol+1,xld+inz+1,n-i-1);;
      inz++;
    }
  }
}

/* ************************************************************
   PROCEDURE utmulx - Compute y = triu(U'*XU), U upper triag, in
   1*n + 2*(n-1) + ... + n*1 mults. Only triu(XU) will be used.
   ************************************************************ */
void utmulx(double *y, const double *u, const double *xu, const mwIndex n)
{
  mwIndex inz,i,j, m;
  const double *ui;

  inz = 0;
  for(j = 0, m = n; j < n; j++, xu += n){
    for(i = 0, ui = u; i <= j; i++, ui += n)
      y[inz++] = realdot(ui,xu,i+1);     /* yij = u(0:i,i)' * xu(0:i,j) */
    inz += --m;                          /* skip tril */
  }
}

/* ************************************************************
   PROCEDURE prpiutmulx - Compute y = triu(U'*XU), U upper triag, in
   4*(1*n + 2*(n-1) + ... + n*1) - n(n-1) mults. Only triu(XU) will be used.
   CAUTION: assumes that diag(imag(U)) == 0.
   ************************************************************ */
void prpiutmulx(double *y, double *ypi, const double *u, const double *upi,
                const double *xu, const double *xupi, const mwIndex n)
{
  mwIndex inz,i,j, m, icol;

  inz = 0;
  for(j = 0, m = n; j < n; j++, xu += n, xupi += n){
    for(i = 0, icol = 0; i <= j; i++, icol += n){
/* yij = u(0:i,i)' * xu(0:i,j) where IM(u(i,i))==0*/
      y[inz] = realdot(u+icol,xu,i+1) + realdot(upi + icol, xupi,i);
      ypi[inz] = realdot(u+icol,xupi,i+1) - realdot(upi + icol, xu, i);
      inz++;
    }
    inz += --m;                          /* skip tril */
  }
}

/* ************************************************************
   PROCEDURE realutxu: Computes U'*X*U in only n*(n+1)^2/2
   multiplications.
   INPUT
     u - n x n full matrix, triu(u) = U.
     x - n x n full symmetric matrix: instead of X*U, we actually
       compute X'*U.
     n - order.
   OUTPUT
     y - n x n full matrix; triu(Y) = triu(U'*X*U).
   WORK
     xu - length n^2 working vector, to store triu(X'*U).
   ************************************************************ */
void realutxu(double *y, const double *u, const double *x,
              const mwIndex n, double *xu)
{
  mwIndex inz, i,j, m, icol,jcol;
/* ------------------------------------------------------------
   Compute xu = triu(X'*U):   1^2 + 2^2 + ... + (n-1)^2 + n^2 mults.
   For each j, let uj point to u(0,j); m is
   remaining length of column j, i.e. m=j.
   ------------------------------------------------------------ */
  inz = 0; jcol = 0;
  for(j = 0, m = n; j < n; j++, jcol += n){
    for(i = 0, icol = 0; i <= j; i++, icol += n)
      xu[inz++] = realdot(x+icol,u+jcol,j+1);  /*xu(i,j) = x(:,i)'*u(0:j,j) */
    inz += --m;                                /* skip tril */
  }
/* ------------------------------------------------------------
   Compute y = triu(U'*XU):  1*n + 2*(n-1) + ... + n*1 mults.
   ------------------------------------------------------------ */
  utmulx(y, u, xu, n);
}

/* ************************************************************
   PROCEDURE prpiutxu: Computes U'*X*U with complex data.
   INPUT
     u - n x n full matrix, triu(u) = U.
     x - n x n full symmetric matrix: instead of X*U, we actually
       compute X'*U.
     n - order.
   OUTPUT
     y - n x n full matrix; triu(Y) = triu(U'*X*U).
   WORK
     xu - length n^2 working vector, to store triu(X'*U).
   ************************************************************ */
void prpiutxu(double *y, double *ypi, const double *u, const double *upi,
              const double *x, const double *xpi,
              const mwIndex n, double *xu)
{
  mwIndex inz, i,j, m, icol,jcol;
  double *xupi;
/* ------------------------------------------------------------
   Partition working array xu into real and imaginary part
   ------------------------------------------------------------ */
  xupi = xu + SQR(n);
/* ------------------------------------------------------------
   Compute xu = triu(X'*U)   ( =triu(XU), since X is Hermitian.)
   For each j, let uj point to u(0,j); m is
   remaining length of column j, i.e. m=j.
   ------------------------------------------------------------ */
  inz = 0; jcol = 0;
  for(j = 0, m = n; j < n; j++, jcol += n){
    for(i = 0, icol = 0; i <= j; i++, icol += n){
/*xu(i,j) = x(:,i)'*u(0:j,j), assume IM u(j,j) == 0 */
      xu[inz] = realdot(x+icol,u+jcol,j+1) + realdot(xpi+icol,upi+jcol,j);
      xupi[inz] = realdot(x+icol,upi+jcol,j) - realdot(xpi+icol,u+jcol,j+1);
      inz++;
    }
    inz += --m;                                /* skip tril */
  }
/* ------------------------------------------------------------
   Compute y = triu(U'*XU)
   ------------------------------------------------------------ */
  prpiutmulx(y, ypi, u, upi, xu, xupi, n);
}

/* ************************************************************
   PROCEDURE partfwsolve -- Solve y from U(1:n,1:n)'*y(1:n) = x(1:n)
      where U is upper-triangular, and of order m >= n.
   INPUT
     u - m x m full matrix; only its upper-triangular entries get used
     x - length n vector.
     m - order of u, m >= n
     n - order of x
   OUTPUT
     y - length n vector, y = (U'\x)_{1:n}.
   ************************************************************ */
void partfwsolve(double *y, const double *u, const double *x,
                const mwIndex m, const mwIndex n)
{
  mwIndex k;
/* ------------------------------------------------------------
   The first equation, u(:,1)'*y=x(1), yields y(1) = x(1)/u(1,1).
   For k = 2:n, we solve
   u(1:k-1,k)'*y(1:k-1) + u(k,k)*y(k) = x(k).
   ------------------------------------------------------------ */
  y[0] = x[0] / u[0];
  u += m;                             /* done with the first column of u */
  for(k = 1; k < n; k++, u += m)
    y[k] = (x[k] - realdot(y,u,k)) / u[k];
}

/* Complex case. Assume IM diag(u) is all-0.
 Solves U' * y = x */
void prpipartfwsolve(double *y,double *ypi, const double *u,const double *upi,
                     const double *x,const double *xpi,
                     const mwIndex m, const mwIndex n)
{
  mwIndex k;
/* ------------------------------------------------------------
   The first equation, u(:,1)'*y=x(1), yields y(1) = x(1)/u(1,1).
   For k = 2:n, we solve
   u(1:k-1,k)'*y(1:k-1) + u(k,k)*y(k) = x(k).
   Assume that IM u(k,k) == 0.
   ------------------------------------------------------------ */
  y[0] = x[0] / u[0];
  ypi[0] = xpi[0] / u[0];
  u += m;                             /* done with the first column of u */
  upi += m;
  for(k = 1; k < n; k++, u += m, upi+=m){
/* y(k) = {x(k) - u(1:k-1,k)'*y(1:k-1)} / ukk */
    y[k] = (x[k] - realdot(u,y,k) - realdot(upi,ypi,k)) / u[k];
    ypi[k] = (xpi[k] - realdot(u,ypi,k) + realdot(upi,y,k)) / u[k];
  }
}

/* ************************************************************
   PROCEDURE partbwsolve -- Solve y from L(j:m,j:m)'*y(j:m) = x(j:m)
      where L is lower-triangular, and of order m > j.
   INPUT
     l - m x m full matrix; only its lower-triangular entries get used
     x - length n vector.
     m - order of l, m > j
     j - last entry of x to be computed (going backwards). If j = 0, then
       all m entries are computed, otherwise only the m-j bottom ones.
       0 <= j < m.
   OUTPUT
     y - length m vector, y[j:m-1] = (L'\x)_{[j:m-1]}.
   ************************************************************ */
void partbwsolve(double *y, const double *l, const double *x,
                 const mwIndex m, const mwIndex j)
{
  mwIndex k;
  const double *lk;
  double ldoty;
/* ------------------------------------------------------------
   The last equation, l(:,m)'*y=x(m), yields y(m) = x(m)/l(m,m).
   For k = m-1:j, we solve
   l(k+1:m,k)'*y(k+1:m) + l(k,k)*y(k) = x(k).
   ------------------------------------------------------------ */
  k = m-1;
  lk = l + m*k;
  ldoty = 0.0;
  for(k = m-1; k > j; k--){
    y[k] = (x[k] - ldoty) / lk[k];
    lk -= m;
    ldoty = realdot(y+k,lk+k,m-k);       /* ldoty = l(k:m,k-1)'*y(k:m) */
  }
  y[j] = (x[j] - ldoty) / lk[j];
}

/* Complex case. Assume IM diag(l) is all-0.
 Solves L' * y = x */
void prpipartbwsolve(double *y,double *ypi, const double *l,const double *lpi,
                     const double *x,const double *xpi,
                     const mwIndex m, const mwIndex j)
{
  mwIndex k;
  const double *lk, *lkpi;
  double ldoty, ldotyim;
/* ------------------------------------------------------------
   The last equation, l(:,m)'*y=x(m), yields y(m) = x(m)/l(m,m).
   For k = m-1:j, we solve
   l(k+1:m,k)'*y(k+1:m) + l(k,k)*y(k) = x(k). Assume l(k,k) == 0.
   ------------------------------------------------------------ */
  k = m-1;
  lk = l + m*k; lkpi = lpi + m*k;
  ldoty = 0.0; ldotyim = 0.0;
  for(k = m-1; k > j; k--){
/* y(k) = (x(k) - l(k+1:m,k)'*y(k+1:m)) / l(k,k) */
    y[k] = (x[k] - ldoty) / lk[k];
    ypi[k] = (xpi[k] - ldotyim) / lk[k];
/* ldoty = l(k:m,k-1)'*y(k:m) */
    lk -= m; lkpi -= m;
    ldoty = realdot(lk+k,y+k,m-k) + realdot(lkpi+k,ypi+k,m-k);;
    ldotyim = realdot(lk+k,ypi+k,m-k) - realdot(lkpi+k,y+k,m-k);;
  }
  y[j] = (x[j] - ldoty) / lk[j];
  ypi[j] = (xpi[j] - ldotyim) / lk[j];
}

/* ************************************************************
   PROCEDURE invutxu: Computes U'\(X/U) in only n*(n+1)^2/2
   multiplications.
   INPUT
     u - m x m full matrix, triu(u) = U.
     x - m x m full symmetric matrix: instead of triu(X/U), we actually
       compute tril(U'\X)' literally.
     m - order.
   OUTPUT
     y - m x m full matrix; triu(Y) = triu(U'\X/U).
   WORK
     xu - length m^2 working vector, to store triu(X/U).
   ************************************************************ */
void invutxu(double *y, const double *u, const double *x,
              const mwIndex m, double *xu)
{
  mwIndex j, jcol;
/* ------------------------------------------------------------
   Compute xu = (X/U)' = U'\X because X is symmetric.
   ------------------------------------------------------------ */
  jcol = 0;
  for(j = 0; j < m; j++, jcol += m)        /* Computes complete U'\X */
    partfwsolve(xu+jcol, u, x+jcol, m,m);
/* ------------------------------------------------------------
   Let triu(xu) = tril(ux)' = triu(X/U)    NB:discards triu(U'\X)
   ------------------------------------------------------------ */
  tril2sym(xu,m);
/* ------------------------------------------------------------
   Compute y = triu(U'\XU):  1*n + 2*(n-1) + ... + n*1 mults.
   viz. fwsolve to j-th entry in column j.
   ------------------------------------------------------------ */
  for(j = 1, jcol = 0; j <= m; j++, jcol += m)
    partfwsolve(y+jcol, u, xu+jcol, m,j);
}

/* ************************************************************
   PROCEDURE prpiinvutxu: Computes U'\(X/U) with complex data.
   INPUT
     u - m x m full matrix, triu(u) = U.
     x - m x m full symmetric matrix: instead of triu(X/U), we actually
       compute tril(U'\X)' literally.
     m - order.
   OUTPUT
     y - m x m full matrix; triu(Y) = triu(U'\X/U).
   WORK
     xu - length 2 * m^2 working vector, to store triu(X/U).
   ************************************************************ */
void prpiinvutxu(double *y,double *ypi, const double *u,const double *upi,
             const double *x,const double *xpi, const mwIndex m, double *xu)
{
  mwIndex j, jcol;
  double *xupi;
/* ------------------------------------------------------------
   Partition xu in real and imaginary part
   ------------------------------------------------------------ */
  xupi = xu + SQR(m);
/* ------------------------------------------------------------
   Compute xu = (X/U)' = U'\X because X is symmetric.
   ------------------------------------------------------------ */
  jcol = 0;
  for(j = 0; j < m; j++, jcol += m)        /* Computes complete U'\X */
    prpipartfwsolve(xu+jcol,xupi+jcol, u,upi, x+jcol,xpi+jcol, m,m);
/* ------------------------------------------------------------
   Let triu(xu) = tril(ux)' = triu(X/U)    NB:discards triu(U'\X)
   ------------------------------------------------------------ */
  tril2herm(xu,xupi,m);
/* ------------------------------------------------------------
   Compute y = triu(U'\XU):  1*n + 2*(n-1) + ... + n*1 mults.
   viz. fwsolve to j-th entry in column j.
   ------------------------------------------------------------ */
  for(j = 1, jcol = 0; j <= m; j++, jcol += m)
    prpipartfwsolve(y+jcol,ypi+jcol, u,upi, xu+jcol,xupi+jcol, m,j);
}

/* ************************************************************
   PROCEDURE invltxl: Computes L'\(X/L) in only n*(n+1)^2/2
   multiplications.
   INPUT
     l - m x m full matrix, tril(l) = L.
     x - m x m full symmetric matrix: instead of tril(X/L), we actually
       compute triu(L'\X)' literally.
     m - order.
   OUTPUT
     y - m x m full matrix; tril(Y) = tril(L'\X/L).
   WORK
     xl - length m^2 working vector, to store tril(X/L).
   ************************************************************ */
void invltxl(double *y, const double *l, const double *x,
              const mwIndex m, double *xl)
{
  mwIndex j, jcol;
/* ------------------------------------------------------------
   Compute xl = tril(X'/L)' = triu(L'\X)
   ------------------------------------------------------------ */
  jcol = 0;
  for(j = 0; j < m; j++, jcol += m)        /* Computes complete L'\X */
    partbwsolve(xl+jcol, l, x+jcol, m,(mwIndex)0);
/* ------------------------------------------------------------
   Let tril(xl) = triu(lx)' = tril(X/L)    NB:discards tril(L'\X)
   ------------------------------------------------------------ */
  triu2sym(xl,m);
/* ------------------------------------------------------------
   Compute y = tril(L'\XL):  1*n + 2*(n-1) + ... + n*1 mults.
   viz. bwsolve to j-th entry in column j=0:m-1.
   ------------------------------------------------------------ */
  for(j = 0, jcol = 0; j < m; j++, jcol += m)
    partbwsolve(y+jcol, l, xl+jcol, m,j);
}

/* complex case. xl is 2*m^2. assume IM diag(l) is all-0. */
void prpiinvltxl(double *y,double *ypi, const double *l,const double *lpi,
                 const double *x,const double *xpi, const mwIndex m, double *xl)
{
  mwIndex j, jcol;
  double *xlpi;
/* ------------------------------------------------------------
   Partition xl in real and imaginary part
   ------------------------------------------------------------ */
  xlpi = xl + SQR(m);
/* ------------------------------------------------------------
   Compute xl = tril(X'/L)' = triu(L'\X)
   ------------------------------------------------------------ */
  jcol = 0;
  for(j = 0; j < m; j++, jcol += m)        /* Computes complete L'\X */
    prpipartbwsolve(xl+jcol,xlpi+jcol, l,lpi, x+jcol,xpi+jcol, m,(mwIndex)0);
/* ------------------------------------------------------------
   Let tril(xl) = triu(lx)' = tril(X/L)    NB:discards tril(L'\X)
   ------------------------------------------------------------ */
  triu2herm(xl,xlpi,m);
/* ------------------------------------------------------------
   Compute y = tril(L'\XL)
   by bwsolve to j-th entry in column j=0:m-1.
   ------------------------------------------------------------ */
  for(j = 0, jcol = 0; j < m; j++, jcol += m)
    prpipartbwsolve(y+jcol,ypi+jcol, l,lpi, xl+jcol,xlpi+jcol, m,j);
}

/* ************************************************************
   PROCEDURE psdscaleK - Computes y = D(d)x over PSD blocks.
   Uses D=U'*U factorization, (transp == 0) Y = UXU'
   or (transp == 1) Y = U'XU.
   INPUT
     x - length lenud input vector.
     ud - Cholesky factor of d for PSD part (after PERM ordering).
     perm - ordering: UD=chol(d(perm,perm)), for numerical stability.
        If perm==NULL, then no reordering is applied.
     cK   - structure describing symmetric cone K.
   OUTPUT
     y - length lenud output vector, y=D(d)x for PSD blocks.
   WORK
     fwork - fwork(2*max(K.s)^2): length 2 * max(rmaxn^2,2*hmaxn^2)
       working vector.
   REMARK lenud := cK.rDim + cK.hDim
   ************************************************************ */
void psdscaleK(double *y, const double *ud, const mwIndex *perm, const double *x,
            const coneK cK, bool transp, double *fwork)
{
  mwIndex k,nk,nksqr;
  double *z, *zpi;
  char use_pivot;
/* ------------------------------------------------------------
   Partition fwork into fwork(psdblk) and z(psdblk), where
   psdblk = max(rmaxn^2,2*hmaxn^2). Let zpi = z+hmaxn^2.
   ------------------------------------------------------------ */
  use_pivot = (perm != (const mwIndex *) NULL);
  z = fwork + MAX(SQR(cK.rMaxn),2*SQR(cK.hMaxn));
  zpi = z + SQR(cK.hMaxn);
/* ------------------------------------------------------------
   PSD: (I) full and !transp
   Y = Ld' * X * Ld. Let Y=X(p,p), where Ld = Ud' (stored in tril(Ud)).
   tril(Y_new) = tril(Ld'* tril(Y*Ld)).
   ------------------------------------------------------------ */
  if(!transp){
    if(use_pivot){            /* with pivoting */
      for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
        nk = cK.sdpNL[k];
        matperm(z,x,perm,nk);
        realltxl(y,ud,z,nk,fwork);
        tril2sym(y,nk);
        nksqr = SQR(nk);
        y += nksqr; ud += nksqr;
        x += nksqr; perm += nk;
      }
      for(; k < cK.sdpN; k++){                     /* complex Hermitian */
        nk = cK.sdpNL[k];
        nksqr = SQR(nk);
        matperm(z,x,perm,nk);
        matperm(zpi,x+nksqr,perm,nk);
        prpiltxl(y,y+nksqr,ud,ud+nksqr,z,zpi,nk,fwork);
        tril2herm(y,y+nksqr,nk);
        nksqr += nksqr;                              /* 2*n^2 for real+imag */
        y += nksqr; ud += nksqr;
        x += nksqr; perm += nk;
      }
    }
    else{ /* without pivoting */
      for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
        nk = cK.sdpNL[k];
        realltxl(y,ud,x,nk,fwork);
        tril2sym(y,nk);
        nksqr = SQR(nk);
        y += nksqr; ud += nksqr;
        x += nksqr;
      }
      for(; k < cK.sdpN; k++){                     /* complex Hermitian */
        nk = cK.sdpNL[k];
        nksqr = SQR(nk);
        prpiltxl(y,y+nksqr,ud,ud+nksqr,x,x+nksqr,nk,fwork);
        tril2herm(y,y+nksqr,nk);
        nksqr += nksqr;                              /* 2*n^2 for real+imag */
        y += nksqr; ud += nksqr;
        x += nksqr;
      }
    }
  }
  else{
/* ------------------------------------------------------------
   (II) transp == 1 then Y = Ud' * X * Ud
   ------------------------------------------------------------ */
    if(use_pivot){            /* with pivoting */
      for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
        nk = cK.sdpNL[k];
        realutxu(z,ud,x,nk,fwork);
        triu2sym(z,nk);
        invmatperm(y,z,perm,nk);            /* Y(perm,perm) = Z */
        nksqr = SQR(nk);
        y += nksqr; ud += nksqr;
        x += nksqr; perm += nk;
      }
      for(; k < cK.sdpN; k++){                     /* complex Hermitian */
        nk = cK.sdpNL[k];
        nksqr = SQR(nk);
        prpiutxu(z,zpi, ud,ud+nksqr,x,x+nksqr,nk,fwork);
        triu2herm(z,zpi,nk);
        invmatperm(y,z,perm,nk);                     /* Y(perm,perm) = Z */
        invmatperm(y+nksqr,zpi,perm,nk);            /* imaginary part */
        nksqr += nksqr;                              /* 2*n^2 for real+imag */
        y += nksqr; ud += nksqr;
        x += nksqr; perm += nk;
      }
    }
    else{ /* without pivoting */
      for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
        nk = cK.sdpNL[k];
        realutxu(y,ud,x,nk,fwork);
        triu2sym(y,nk);
        nksqr = SQR(nk);
        y += nksqr; ud += nksqr;
        x += nksqr;
      }
      for(; k < cK.sdpN; k++){                     /* complex Hermitian */
        nk = cK.sdpNL[k];
        nksqr = SQR(nk);
        prpiutxu(y,y+nksqr, ud,ud+nksqr,x,x+nksqr,nk,fwork);
        triu2herm(y,y+nksqr,nk);
        nksqr += nksqr;                              /* 2*n^2 for real+imag */
        y += nksqr; ud += nksqr;
        x += nksqr;
      }
    }
  }
}

#ifdef SEDUMI_OLD
/* ************************************************************
   PROCEDURE scaleK - Computes y = D(d)x.
   For PSD, uses D=U'*U factorization, (transp == 0) Y = UXU'
   or (transp == 1) Y = U'XU.
   INPUT
     x - length N(K) input vector.
     d - scaling vector, only LP and Lorentz part needed.
     ud - Cholesky factor of d for PSD part (after PERM ordering).
     qdetd - sqrt(det(d)) for Lorentz part.
     perm - ordering: UD=chol(d(perm,perm)), for numerical stability.
     cK   - structure describing symmetric cone K.
     invdx - length cK.qDim vector containing D(d)\x for Lorentz part.
       This is optional. If invdx == NULL, then not used.
   OUTPUT
     y - length N(K) output vector, y=D(d)x.
     dmult - If !NULL then lorN-vector containing mu[k], such that
       (D(d)x)_k = y_k + mu[k] * d_k, where "_k" is the kth Lorentz block.
   WORK
     fwork - fwork(2*max(K.s)^2): length 2 * max(rmaxn^2,2*hmaxn^2)
       working vector.
   ************************************************************ */
void scaleK(double *y, double *dmult, const double *d, const double *ud,
            const double *qdetd, const mwIndex *perm, const double *x,
            const coneK cK, const char transp, const double *invdx,
            double *fwork)
{
  mwIndex k,nk;
  double detdk;
/* ------------------------------------------------------------
   LP: y = d .* x
   ------------------------------------------------------------ */
  realHadamard(y, d,x,cK.lpN);
  y += cK.lpN;             /* Next, point to lorentz & sdp blocks */
  d += cK.lpN; x += cK.lpN;
/* ------------------------------------------------------------
   LORENTZ (1/3): y = D(d) x
   ------------------------------------------------------------ */
  if(dmult == (double *) NULL)
    for(k = 0; k < cK.lorN; k++){
      nk = cK.lorNL[k];
      qlmul(y,d,x,qdetd[k],nk);
      y += nk; d += nk; x +=nk;
    }
  else
/* ------------------------------------------------------------
   LORENTZ (2/3): D(d) x = y + dmult * d. This storage scheme avoids
   cancelation in y.
   ------------------------------------------------------------ */
    if(invdx == (const double *) NULL)
      for(k = 0; k < cK.lorN; k++){
        nk = cK.lorNL[k];
        dmult[k] = qscale(y,d,x,qdetd[k],nk);
        y += nk; d += nk; x +=nk;
      }
    else{
/* ------------------------------------------------------------
   LORENTZ (3/3): D(d) x = D(d^2) invdx = y + dmult * d
   USES D(d^2)invdx = (d'*invdx) * d + det(d) * [-invdx(1); invdx(2:nk)]
   We let y = det(d) * [-invdx(1); invdx(2:nk)].
   ------------------------------------------------------------ */
      for(k = 0; k < cK.lorN; k++){
        nk = cK.lorNL[k];
        dmult[k] = realdot(d,invdx,nk);
        detdk = SQR(qdetd[k]);
        y[0] = - detdk * invdx[0];
        scalarmul(y+1, detdk,invdx+1,nk-1);
        y += nk; d += nk; invdx +=nk;
      }
      x += cK.qDim;         /* point beyond Lorentz */
    }
/* ------------------------------------------------------------
   PSD scale
   ------------------------------------------------------------ */
  psdscaleK(y, ud, perm, x, cK, transp, fwork);
}
#endif
