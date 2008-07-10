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
#include "blksdp.h"
#include "reflect.h"

/* ************************************************************
   PROCEDURE elqxq - Compute (I+c*c'/beta) * [X1, X2 * (I+c*c'/beta)];
     only tril is computed, and only tril(X2) is used.
   INPUT
     beta - Householder scalar coefficient
     c    - length m-n1 elementary Householder vector
     m    - length of columns in X. m >= n1.
     n1   - number of X1-columns; the order of Householder reflection will
         be m-n1.
   UPDATED
     x -  (m x m)-n1. On output, Xnew = (I+c*c'/beta)* [X1, X2*(I+c*c'/beta)]
       only tril is computed. We treat x as a (m-n1)*m matrix, but we need to add
       m (instead of m-n1) to get to the next column.
   WORK
     y - length m working vector, for storing c'*[X1, X2].
   ************************************************************ */
void elqxq(double *x, const double beta, const double *c,
           const mwIndex n1, const mwIndex m, double *y)
{
  mwIndex j,n2;
  double *xj, *y2, *x2;
  double alpha;
  n2 = m - n1;          /* order of x2 */
  y2 = y + n1;

/* ------------------------------------------------------------
   Compute y1 = c'*X1
   ------------------------------------------------------------ */
  for(j = 0, xj = x; j < n1; j++, xj += m)
    y[j] = realdot(c, xj, n2);
  x2 = xj;


/* ------------------------------------------------------------
   Compute y2 = c'*tril(X2); y2 += tril(X2,-1)*c, SO THAT y2 = c'*X2SYM.
   ------------------------------------------------------------ */
  for(j = 0; j < n2; j++, xj += m+1)
    y2[j] = realdot(c+j, xj, n2-j);

  for(j = 1, xj = x2+1; j < n2; j++, xj += m+1)
    addscalarmul(y2+j, c[j-1], xj, n2-j);
  

/* ------------------------------------------------------------
   Below-diag block: let X1 += c*y1' / beta
   ------------------------------------------------------------ */
  for(j = 0, xj = x; j < n1; j++, xj += m)
    addscalarmul(xj, y[j] / beta, c, n2);
  

/* ------------------------------------------------------------
   Lower-Right block: X2 += c*((y2'*c/beta) * c' + y2')/beta + y2*c'/beta
   ------------------------------------------------------------ */
  alpha = realdot(y2, c, n2) / beta;
 
  for(j = 0, xj = x2; j < n2; j++, xj += m+1){
    addscalarmul(xj, (alpha * c[j] + y2[j])/beta, c+j, n2-j);
    addscalarmul(xj, c[j] / beta, y2+j, n2-j);
  }

}

/* ************************************************************
   PROCEDURE prpielqxq - Compute (I+c*c'/beta) * [X1, X2 * (I+c*c'/beta)];
     only tril is computed, and only tril(X2) is used.
   INPUT
     beta - Householder scalar coefficient (real length m)
     c,cpi - length m-n1 elementary Householder vector
     m    - length of columns in X. m >= n1.
     n1   - number of X1-columns; the order of Householder reflection will
         be m-n1.
   UPDATED
     x,xpi - 2*((m x m)-n1). On output,
       Xnew = (I+c*c'/beta)* [X1, X2*(I+c*c'/beta)]
       only tril is computed.  We treat x as a (m-n1)*m matrix, but we need to add
       m (instead of m-n1) to get to the next column.
   WORK
     y - length 2*m working vector, for storing c'*[X1, X2].
   ************************************************************ */
void prpielqxq(double *x, double *xpi, const double beta, const double *c,
               const double *cpi, const mwIndex n1, const mwIndex m, double *y)
{
  mwIndex j,n2;
  double *xj,*xjpi, *y2, *x2, *x2pi, *ypi, *y2pi;
  double alpha;
  n2 = m - n1;          /* order of x2 */
/* ------------------------------------------------------------
   Partition y into y(n1), y2(n2); ypi(n1), ypi(n2).
   ------------------------------------------------------------ */
  ypi = y + m;
  y2 = y + n1;
  y2pi = ypi + n1;
/* ------------------------------------------------------------
   Compute y1 = c'*X1.   NOTE: y1 is 1 x n1 row-vector.
   ------------------------------------------------------------ */
  for(j = 0, xj = x, xjpi = xpi; j < n1; j++, xj += m, xjpi += m){
    y[j] = realdot(c, xj, n2) + realdot(cpi, xjpi, n2);
    ypi[j] = realdot(c, xjpi, n2) - realdot(cpi, xj, n2);
  }
  x2 = xj; x2pi = xjpi;
/* ------------------------------------------------------------
   Compute y2 = c'*tril(X2)
   ------------------------------------------------------------ */
  for(j = 0; j < n2; j++, xj += m+1, xjpi += m+1){
/* y2(j) = c(j:n2)'*x(j:n2,j) */
    y2[j] = realdot(c+j, xj, n2-j) + realdot(cpi+j, xjpi, n2-j);
    y2pi[j] = realdot(c+j, xjpi, n2-j) - realdot(cpi+j, xj, n2-j);
  }
/* ------------------------------------------------------------
   Let y2 += (tril(X2,-1)*c)' = c'*triu(X2,1),
   SO THAT y2 = c'*X2, with X2 symmetric. NOTE: y2 = 1 x n2 row-vector.
   ------------------------------------------------------------ */
  for(j = 1, xj = x2+1, xjpi = x2pi+1; j < n2; j++, xj += m+1, xjpi += m+1){
/* y2(j+1:n2) += conj(x(j+1:n2,j) * c(j)) */
    addscalarmul(y2+j, c[j-1], xj, n2-j);        /* RE */
    addscalarmul(y2+j, -cpi[j-1], xjpi, n2-j);
    addscalarmul(y2pi+j, -cpi[j-1], xj, n2-j);   /* -IM, i.e. conj */
    addscalarmul(y2pi+j, -c[j-1], xjpi, n2-j);
  }
/* ------------------------------------------------------------
   Below-diag block: let X1 += c*y1 / beta, where y1 = c'*X1.
   NOTE: y1 is 1 x n1 row-vector.
   This completes X1_new = (I+c*c'/beta) * X1_old.
   ------------------------------------------------------------ */
  for(j = 0, xj = x, xjpi = xpi; j < n1; j++, xj += m, xjpi += m){
/* x(:,j) += c * y1(j) / beta */
    addscalarmul(xj, y[j] / beta, c, n2);          /* RE */
    addscalarmul(xj, -ypi[j] / beta, cpi, n2);
    addscalarmul(xjpi, y[j] / beta, cpi, n2);          /* IM */
    addscalarmul(xjpi, ypi[j] / beta, c, n2);
  }
/* ------------------------------------------------------------
   Lower-Right block: X2 += c*((y2*c/beta) * c' + y2)/beta + y2'*c'/beta
   where y2 = c'*X2.
   This completes X2new = (I+c*c'/beta) * X2 * (I+c*c'/beta)
   NOTE: since X2 = X2', we have y2*c = c'*X2*c is real.
   ------------------------------------------------------------ */
/* alpha = y2 * c / beta, which is real.*/
  alpha = (realdot(y2, c, n2) - realdot(y2pi, cpi, n2)) / beta;
  for(j = 0, xj = x2, xjpi = x2pi; j < n2; j++, xj += m+1, xjpi += m+1){
/* x2(j:n2,j) += c(j:n2) * (alpha*conj(c(j)) + y2(j))/beta */
    addscalarmul(xj, (alpha * c[j] + y2[j])/beta, c+j, n2-j);
    addscalarmul(xj, (alpha * cpi[j] - y2pi[j])/beta, cpi+j, n2-j);
    addscalarmul(xjpi, (alpha * c[j] + y2[j])/beta, cpi+j, n2-j);
    addscalarmul(xjpi, (y2pi[j] - alpha * cpi[j])/beta, c+j, n2-j);
/* x2(j:n2,j) += conj(c(j)*y2(j:n2)) / beta */
    addscalarmul(xj, c[j] / beta, y2+j, n2-j);            /* RE */
    addscalarmul(xj, -cpi[j] / beta, y2pi+j, n2-j);
    addscalarmul(xjpi, -c[j] / beta, y2pi+j, n2-j);       /* -IM, i.e. conj */
    addscalarmul(xjpi, -cpi[j] / beta, y2+j, n2-j);
  }
}

/* ************************************************************
   PROCEDURE qtxq - computes tril(Qb' * X * Qb)
    Here, Qb = Q_1*Q_2*..*Q_{m-1}, where each Q_i is a Householder reflection.
    (Qb is from a Qb * R decomposition.)
   INPUT
     beta - length m vector
     c    - m x m matrix, lower triangular gives Householder reflections
     m    - order
   UPDATED
     x -  m x m. On output, Xnew = Qb' * X * Qb
       This means: start with order m reflection, up to order 2 reflection.
   WORK
     fwork - length m working vector.
   ************************************************************ */
void qtxq(double *x, const double *beta, const double *c,
          const mwIndex m, double *fwork)
{
  mwIndex k, inz;

  inz = 0;
/* ------------------------------------------------------------
   For each k, c[inz] = c(k,k), the top of the lower-right block,
   x[k] is start of k-th row in k.
   ------------------------------------------------------------ */
  for(k = 0; k < m-1; k++, inz += m+1)
    elqxq(x + k, -beta[k], c + inz, k, m, fwork);
}

/* ************************************************************
   PROCEDURE prpiqtxq - computes tril(Qb' * X * Qb)
    Here, Qb = Q_1*Q_2*..*Q_{m-1}*diag(q(:,m)), where each Q_i is a
    Householder reflection, and q(:,m) is a complex sign-vector.
    (Qb is from a Qb * R decomposition.)
   INPUT
     beta - length m vector (real)
     c,cpi - m x m matrix, lower triangular gives Householder reflections
     m    - order
   UPDATED
     x,xpi -  m x m. On output, Xnew = Qb' * X * Qb
       This means: start with order m reflection, up to order 2 reflection.
   WORK
     fwork - length 2*m working vector.
   ************************************************************ */
void prpiqtxq(double *x, double *xpi, const double *beta, const double *c,
              const double *cpi, const mwIndex m, double *fwork)
{
  mwIndex i,k, inz;
  const double *qsgn, *qsgnpi;
  double qk, qkim, qij,qijim, xij;

/* ------------------------------------------------------------
   START by Householder transformations:
   For each k, c[inz] = c(k,k), the top of the lower-right block,
   x[k] is start of k-th row in k.
   ------------------------------------------------------------ */
  inz = 0;
  for(k = 0; k < m-1; k++, inz += m+1)
    prpielqxq(x + k,xpi + k, -beta[k], c + inz,cpi + inz, k, m, fwork);
/* ------------------------------------------------------------
   FINISH BY COMPLEX SIGNING.
   Let qsgn = c(:,m) be the sign vector. Then let
   X_new = diag(qsgn)' * X * diag(qsgn) = conj(qsign(i)) * qsign(j) * x_ij
   for i > j. The diagonal is not affected, since |qsign(i)| = 1.
   ------------------------------------------------------------ */
  qsgn = c + m * (m-1);
  qsgnpi = cpi + m * (m-1);
  inz = 0;
  for(k = 0; k < m-1; k++){
    qk = qsgn[k]; qkim = qsgnpi[k];
    inz += k+1;                       /* point below diagonal */
    for(i = k+1; i < m; i++){
/* qij = conj(qsign(i)) * qsign(k) */
      qij = qsgn[i]*qk + qsgnpi[i] * qkim;
      qijim = qsgn[i]*qkim - qsgnpi[i] * qk;
/* xij *= qij */
      xij = x[inz] * qij - xpi[inz] * qijim;
      xpi[inz] = x[inz] * qijim + xpi[inz] * qij;
      x[inz] = xij;
      inz++;
    }
  }
}

/* ************************************************************
   PROCEDURE qxqt - computes tril(Qb * X * Qb')
    Here, Qb = Q_1*Q_2*..*Q_{m-1}, where each Q_i is a Householder reflection.
    (Qb is from a Qb * R decomposition.)
   INPUT
     beta - length m vector
     c    - m x m matrix, lower triangular gives Householder reflections
     m    - order
   UPDATED
     x -  m x m. On output, Xnew = Qb * X * Qb'
       This means: start with order 2 reflection, up to order m reflection.
   WORK
     fwork - length m working vector.
   ************************************************************ */
void qxqt(double *x, const double *beta, const double *c,
          const mwIndex m, double *fwork)
{
  mwIndex k, inz;
    mxAssert(m>1,"");
  inz = SQR(m) - (m+2);

/* ------------------------------------------------------------
   For each k, c[inz] = c(k,k), the top of the lower-right block,
   x[k] is start of k-th row in k.
   ------------------------------------------------------------ */
  for(k = m-1; k > 0; k--, inz -= m+1){
    elqxq(x + k-1, -beta[k-1], c + inz, k-1, m, fwork);
  }
  /*elqxq(x , -beta[0], c , 0, m, fwork);*/

}

/* ************************************************************
   PROCEDURE prpiqxqt - computes tril(Qb * X * Qb')
    Here, Qb = Q_1*Q_2*..*Q_{m-1}*diag(q(:,m)), where each Q_i is a
    Householder reflection, and q(:,m) is a complex sign-vector.
    (Qb is from a Qb * R decomposition.)
   INPUT
     beta - length m vector (real)
     c,cpi - m x m matrix, lower triangular gives Householder reflections
     m    - order
   UPDATED
     x,xpi -  m x m. On output, Xnew = Qb * X * Qb'
       This means: start with order 2 reflection, up to order m reflection.
   WORK
     fwork - length 2*m working vector.
   ************************************************************ */
void prpiqxqt(double *x, double *xpi, const double *beta, const double *c,
              const double *cpi, const mwIndex m, double *fwork)
{
  mwIndex i,k, inz;
  const double *qsgn, *qsgnpi;
  double qk, qkim, qij,qijim, xij;

/* ------------------------------------------------------------
   START BY COMPLEX SIGNING.
   Let qsgn = c(:,m) be the sign vector. Then let
   X_new = diag(qsgn) * X * diag(qsgn)' = qsign(i) * conj(qsign(j)) * x_ij
   for i > j. The diagonal is not affected, since |qsign(i)| = 1.
   ------------------------------------------------------------ */
  qsgn = c + m * (m-1);
  qsgnpi = cpi + m * (m-1);
  inz = 0;
  for(k = 0; k < m-1; k++){
    inz += k+1;                       /* point below diagonal */
    qk = qsgn[k]; qkim = qsgnpi[k];
    for(i = k+1; i < m; i++){
/* qij = conj(qsign(i)) * qsign(k) */
      qij = qsgn[i]*qk + qsgnpi[i] * qkim;
      qijim = qsgn[i]*qkim - qsgnpi[i] * qk;
/* xij *= conj(qij) */
      xij = x[inz] * qij + xpi[inz] * qijim;
      xpi[inz] = xpi[inz] * qij - x[inz] * qijim;      /* conj qij */
      x[inz] = xij;                             /* WRITE signed x-values */
      inz++;
    }
  }
/* ------------------------------------------------------------
   FINISH by Householder transformations:
   For each k, c[inz] = c(k,k), the top of the lower-right block,
   x[k] is start of k-th row in k.
   ------------------------------------------------------------ */
  inz = SQR(m) - (m+2);
  for(k = m-1; k > 0; k--, inz -= m+1)
    prpielqxq(x + k-1,xpi + k-1, -beta[k-1], c + inz,cpi+inz, k-1, m, fwork);
}
