/*
  [q,r] = qrK(x,K);

  * x should consist only of SDP part (length sum(K.s.^2)).
  * q is in product form: the n-1 factors are:
    I- q(k:n,k)*q(k:n,k)' / beta(k).
    For real PSD blocks, beta[k] = q(k,n)
    For complex PSD blocks, beta[k] = q(k,2*n+1), and q is n * (2n+1).
  * tril(r) will be all-0.

 qrK: Q(beta)*R factorization

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

#define Q_OUT myplhs[0]
#define R_OUT myplhs[1]

#define NPAROUT 2

#define X_IN prhs[0]
#define K_IN prhs[1]
#define NPARIN 2

/* ************************************************************
   PROCEDURE isconjhadamul - Let y = conj(x).*y
   ************************************************************ */
void isconjhadamul(double *y, double *ypi, const double *x,const double *xpi,
                   const mwIndex n)
{
  mwIndex i;
  double yi;
  for(i = 0; i < n; i++){
    yi = x[i] * y[i] + xpi[i] * ypi[i];
    ypi[i] = x[i] * ypi[i] - xpi[i] * y[i];
    y[i] = yi;
  }
}

/* ************************************************************
   PROCEDURE qrfac - QR factorization for nxn matrix.
   INPUT
     n - order of matrix to be factored
   UPDATED
     u - Full nxn. On input, u is matrix to be factored. On output,
       triu(u) = uppertriangular factor;
       tril(u,-1) = undefined.
   OUTPUT
     beta - length n vector. kth Householder reflection is
        Qk = I-qk*qk' / beta[k],   where qk = q(k:n-1,k).
     q - n x (n-1) matrix; each column is a Householder reflection.
   ************************************************************ */
void qrfac(double *beta, double *q, double *u, const mwIndex n)
{
  mwIndex i,k, kcol, nmink, icol;
  double dk, betak, qkui, qkk;

  for(k = 0, kcol = 0; k < n-1; k++, kcol += n+1){

/* ------------------------------------------------------------
   kth Householder reflection:
   dk = sign(xkk) * ||xk(k:n)||,
   qk(k+1:n) = x(k+1:n); qkk = xkk+dk, betak = dk*qkk, ukk = -dk.
   ------------------------------------------------------------ */
       
    qkk = u[kcol];
    dk = SIGN(qkk) * sqrt(realssqr(u+kcol,n-k));
    memcpy(q + kcol+1, u+kcol+1, (n-k-1) * sizeof(double));
    qkk += dk;
    betak = dk * qkk;
    q[kcol] = qkk;
    if(betak == 0.0)              /* If xk is all-0 then set beta = 1. */
      betak = 1.0;
    beta[k] = betak;
    u[kcol] = -dk;
/* ------------------------------------------------------------
   Reflect columns k+1:n-1, i.e.
   xi -= (qk'*xi / betak) * qk, where xi = x(k:n-1, i).
   ------------------------------------------------------------ */
    nmink = n-k;
    betak = -betak;
  
    for(i = k + 1, icol = kcol + n; i < n; i++, icol += n){
      qkui = realdot(q+kcol, u+icol, nmink);
      addscalarmul(u+icol, qkui/betak, q+kcol, nmink);

    }
  }
}

/* complex Hermitian: */
/* ************************************************************
   PROCEDURE prpiqrfac - QR factorization for nxn matrix.
   INPUT
     n - order of matrix to be factored
   UPDATED
     u - Full nxn. On input, u is matrix to be factored. On output,
       triu(u) = uppertriangular factor;
       tril(u,-1) = undefined.
   OUTPUT
     beta - length n vector. kth Householder reflection is
        Qk = I-qk*qk' / beta[k],   where qk = q(k:n-1,k).
     q,qpi - n x n matrix; each of the first n-1 columns is a Householder
       reflection. The n-th column gives Qn = diag(q(:,n)), which is NOT
       Hermitian (viz. diag complex rotations). We have
       u_IN = Q_1*Q_2*.. *Q_n * triu(u_OUT),
       u_OUT = Q_n' * Q_{n-1}* .. * Q_2*Q_1*u_IN.
   ************************************************************ */
void prpiqrfac(double *beta, double *q, double *qpi, double *u,
               double *upi, const mwIndex n)
{
  mwIndex i,j,k, kcol, nmink, icol;
  double betak, qkui,qkuiim, absxkk, normxk, xkk,xkkim;
  double *ui,*uipi, *qk, *qkpi;

  for(k = 0, kcol = 0; k < n-1; k++, kcol += n+1){
    qk = q+kcol;   
    qkpi = qpi + kcol;
/* ------------------------------------------------------------
   kth Householder reflection:
   Set absxkk = |xkk| and normxk = norm(xk(k:n)), then
   ukk = -sign(xkk) * normxk,     qkk = xkk - ukk,   qk(k+1:n) = xk(k+1:n).
   Remark: sign(xkk) := xkk/|xkk|, a complex number.
   ------------------------------------------------------------ */
    xkk = u[kcol];
    xkkim = upi[kcol];
    absxkk = SQR(xkk) + SQR(xkkim);
    normxk = absxkk + realssqr(u+kcol+1,n-k-1) + realssqr(upi+kcol+1,n-k-1);
    memcpy(qk+1, u+kcol+1, (n-k-1) * sizeof(double));      /* real */
    memcpy(qkpi+1, upi+kcol+1, (n-k-1) * sizeof(double));    /* imag */
    absxkk = sqrt(absxkk);
    normxk = sqrt(normxk);
    if(absxkk > 0.0){
      u[kcol] = -(xkk / absxkk) * normxk;     /* ukk = -sign(xkk) * normxk */
      upi[kcol] = -(xkkim / absxkk) * normxk;
    }
    else
      u[kcol] = -normxk;                      /* sign(0) := 1 */
    qk[0]   = xkk - u[kcol];                  /* qkk = xkk - ukk */
    qkpi[0] = xkkim - upi[kcol];
/* ------------------------------------------------------------
   betak = normxk * (normxk + absxkk)
   EXCEPTION: if xk is all-0 then set beta = 1 (to avoid division by 0).
   ------------------------------------------------------------ */
    if(normxk == 0.0)
      betak = 1.0;
    else
      betak = normxk * (normxk + absxkk);
    beta[k] = betak;
/* ------------------------------------------------------------
   Reflect columns k+1:n-1, i.e.
   xi -= (qk'*xi / betak) * qk, where xi = x(k:n-1, i).
   ------------------------------------------------------------ */
    nmink = n-k;
    for(i = k + 1, icol = kcol + n; i < n; i++, icol += n){
      ui = u + icol; uipi = upi+icol;
      qkui   = realdot(qk, ui, nmink) + realdot(qkpi, uipi, nmink);
      qkuiim = realdot(qk, uipi, nmink) - realdot(qkpi, ui, nmink);
      qkui /= betak;
      qkuiim /= betak;
/* for all j, we have x(j,i) -= (qkui + i * qkuiim) * (qk[j] + i*qkpi[j]) */
      for(j = 0; j < nmink; j++){
        ui[j] -= qkui * qk[j] - qkuiim * qkpi[j];
        uipi[j] -= qkui * qkpi[j] + qkuiim * qk[j];
      }
    }
  } /* k = 0:n-2 */
/* ------------------------------------------------------------
   The Q*R decomposition is now done, but IM diag(u) may be nonzero.
   Therefore, we multiply each row with conj(sign(u_ii)) = conj(u_ii)/|u_ii|.
   Let q(1:n,n) =  sign(diag(u))
   ------------------------------------------------------------ */
  mxAssert(n>=0,"");
  if(n > 0){
    kcol = n * (n-1);     /* sign column q(:,n) */
    qk = q + kcol;
    qkpi = qpi + kcol;
/* Let icol point to (i,i) entry */
    for(i = 0, icol = 0; i < n; i++, icol += n+1){
      xkk = u[icol];                            /* get u(i,i) */
      xkkim = upi[icol];
      absxkk = sqrt(SQR(xkk) + SQR(xkkim));
      qk[i] = xkk / absxkk;                     /* q(i) = sign(u(i,i)) */
      qkpi[i] = xkkim / absxkk;
      u[icol] = absxkk;                         /* new u(i,i) = |uOLD(i,i)| */
      upi[icol] = 0.0;
    }
/* ------------------------------------------------------------
   Let Unew = Q_n'*Uold, i.e. u(i,k) = conj(qn(i)) * uOLD(i,k), i < k
   ------------------------------------------------------------ */
    for(k = 1, kcol = n; k < n; k++, kcol += n)
      isconjhadamul(u+kcol, upi+kcol, qk,qkpi,k);
  }
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   [beta,U,d,perm] = qrpfacK(x,K)
   ************************************************************ */
void mexFunction( int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  mxArray *myplhs[NPAROUT];
  coneK cK;
  mwIndex i,k,nk,nksqr, sdpdim, qsize;
  double *q, *r, *betak;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "qrK requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "qrK produces less output arguments");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Compute statistics: sdpdim = rdim+hdim, qsize = sdpdim + hLen.
   ------------------------------------------------------------ */
  sdpdim = cK.rDim + cK.hDim;
  qsize = sdpdim + cK.hLen;
/* ------------------------------------------------------------
   Check input vector x.
   ------------------------------------------------------------ */
  mxAssert(mxGetM(X_IN) * mxGetN(X_IN) == sdpdim, "size mismatch x");
/* ------------------------------------------------------------
   Allocate output Q(qsize), R(sdpdim)
   and let r = x.
   ------------------------------------------------------------ */
  Q_OUT = mxCreateDoubleMatrix(qsize, (mwSize)1, mxREAL);
  q = mxGetPr(Q_OUT);
  R_OUT = mxCreateDoubleMatrix(sdpdim, (mwSize)1, mxREAL);
  r = mxGetPr(R_OUT);
  memcpy(r, mxGetPr(X_IN), sdpdim * sizeof(double));
/* ------------------------------------------------------------
   The actual job is done here:
   ------------------------------------------------------------ */
  for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
    nk = (mwIndex) cK.sdpNL[k];
    nksqr = SQR(nk);
    qrfac(q+nksqr-nk,q,r, nk);
    r += nksqr; 
    q += nksqr;
  }
  for(; k < cK.sdpN; k++){                      /* complex Hermitian */
    nk = (mwIndex) cK.sdpNL[k];
    nksqr = SQR(nk);
    betak = q + 2*nksqr;
    prpiqrfac(betak,q,q+nksqr, r,r+nksqr, nk);
    nksqr += nksqr;
    r += nksqr; 
    q += nksqr + nk;               /* nk for betak */
  }
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  i = MAX(nlhs, 1);
  memcpy(plhs,myplhs, i * sizeof(mxArray *));
  for(; i < NPAROUT; i++)
    mxDestroyArray(myplhs[i]);
}
