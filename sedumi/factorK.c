/*
  [qdetx,ux,ispos,perm] = factorK(x,K);         // with pivoting
  [qdetx,ux,ispos] = factorK(x,K);              // without pivoting

  cholpfac: UX'*UX Cholesky factorization with pivoting for PSD,
  and stable qdet(x) for LORENTZ.
  If ispos = 1 then success
     ispos = 0 then x not in K.

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

#define QDETX_OUT myplhs[0]
#define UX_OUT myplhs[1]
#define ISPOS_OUT myplhs[2]
#define NPARNOPERM 3
#define PERM_OUT myplhs[3]
#define NPAROUT 4

#define X_IN prhs[0]
#define K_IN prhs[1]
#define NPARIN 2


/* ============================================================
   LORENTZ DETERMINANT
   ============================================================ */

/* ************************************************************
   PROCEDURE qdet  -  computes det(x) = lab1 * lab2 w.r.t. Lorentz cone
   INPUT:
     x - full n x 1
     n - length of x
   RETURNS:
     determinant qdet(x).
   ************************************************************ */
double qdet(const double *x,const int n)
{
  double y;
/* ------------------------------------------------------------
   qdet(x) = x'*J*x / 2 = (x_1^2 - |x_2|^2)/2
   For stability, we evaluate it is (x_1+|x_2|)(x_1-|x_2|)/2
   ------------------------------------------------------------ */
  y = sqrt(realssqr(x+1,n-1));
  return ((x[0]+y) * (x[0]-y)) / 2;
}

/* ============================================================
   PSD: CHOLESKY FACTORIZATION
   ============================================================ */

/* ************************************************************
   PROCEDURE cholnopiv - U'*U factorization for nxn matrix,
     without pivotting.
   INPUT
     n - order of matrix to be factored
   UPDATED
     u - Full nxn. Input: Matrix to be factored.
       Output: Cholesky factor, X=triu(U)'*triu(U). 
       NOTE: tril(U,-1) is not affected by this function.
   RETURNS:
     0 = "success", 1 = "X is NOT positive definite"
   ************************************************************ */
int cholnopiv(double *u,const int n)
{
  int i,j;
  double uij,ujj;
  double *ui,*uj;

  if(n < 1)
    return 0;
/* ------------------------------------------------------------
   Solve the columns of U, for j=0:n-1
   ------------------------------------------------------------ */
  for(uj = u, j = 0; j < n; j++, uj+=n){
/* ------------------------------------------------------------
   Solve "uij" from the identity
   uii * uij = xij - u(1:i-1,i)'*u(1:i-1,j)
   ------------------------------------------------------------ */
    for(ui = u, i = 0, ujj = 0.0; i < j; i++, ui+=n){
      uij = (uj[i] = (uj[i] - realdot(ui,uj,i)) / ui[i]);
      ujj += SQR(uij);
    }
    ujj = uj[j] - ujj;
/* ------------------------------------------------------------
   By now, "ujj" should contain the final u(j,j)^2. Check whether
   it is positive. If not, then X was not p.d., thus fail.
   ------------------------------------------------------------ */
    if(ujj <= 0.0)
      return 1;               /* X is not positive definite */
    uj[j] = sqrt(ujj);
  }
  return 0;   /* success */
}

/* ************************************************************
   PROCEDURE prpicholnopiv  -  Computes triu matrix U s.t. U'*U = X.
     U, X in SeDuMi's complex format: X = [RE X, IM X]. No pivoting.
   INPUT:
     x - full 2*(n x n), should be Hermitian.
     n - order of u, x.
   UPDATED:
     u,upi - Both full nxn. Input: Matrix to be factored (real, imaginary).
       Output: Cholesky factor, X=U'*U with U=triu(u) + i*triu(upi,1). 
       NOTE: tril(u,-1) and tril(upi) are not affected by this function.
   RETURNS:
     0 = "success", 1 = "X is NOT positive definite"
   ************************************************************ */
int prpicholnopiv(double *u, double *upi, const int n)
{
  int i,j;
  double uii,uij,ujj;
  double *ui,*uipi, *uj, *ujpi;

  if(n < 1)
    return 0;
/* ------------------------------------------------------------
   Solve the columns of U, for j=0:n-1
   ------------------------------------------------------------ */
  for(uj=u, ujpi=upi, j = 0; j < n; j++, uj+=n, ujpi+=n){
    ujj = 0.0;
    for(ui=u, uipi=upi, i = 0; i < j; i++, ui+=n, uipi+=n){
/* ------------------------------------------------------------
   Solve "uij" from the identity
   uii * uij = xij - u(1:i-1,i)'*u(1:i-1,j)
   ------------------------------------------------------------ */
      uii = ui[i];
      uij = uj[i];                                     /* real part */
      uij -= realdot(ui,uj,i) + realdot(uipi,ujpi,i);
      uj[i] = (uij /= uii);
      ujj += SQR(uij);
      uij = ujpi[i];                                  /* imaginary part */
      uij += realdot(uipi,uj,i) - realdot(ui,ujpi,i);
      ujpi[i] = (uij /= uii);
      ujj += SQR(uij);
    }
    ujj = uj[j] - ujj;
/* ------------------------------------------------------------
   By now, "ujj" should contain the final u(j,j)^2. Check whether
   it is positive. If not, then X was not p.d., thus fail.
   ------------------------------------------------------------ */
    if(ujj <= 0.0)
      return 1;               /* X is not positive definite */
    uj[j] = sqrt(ujj);
  }
  return 0;
}

/* ************************************************************
   PROCEDURE cholpivot - U'*U factorization for nxn matrix,
     with column pivotting. Yields U with 1 >= \|U(i,i+1:n)\| forall i.
   INPUT
     n - order of matrix to be factored
     x - Full nxn. Matrix to be factored.
   OUTPUT
     u - Full nxn, Cholesky factor: X=U'*U with U(:,perm)
       upper triangular.
     perm - Column permutation for maximal stability.
   WORK
     d - length n vector; the positive diagonal diag(U(:,perm)).^2,
         has decreasing order.
   RETURNS:
     0 = "success", 1 = "X is NOT positive definite"
   ************************************************************ */
int cholpivot(double *u,int *perm, const double *x, const int n, double *d)
{
  int i,j,k,imax, icol;
  double *uk, *rowuk;
  double dk, uki;
  const double *xk;
/* ------------------------------------------------------------
   Initialize: d = diag(x), perm = 0:n-1.
   ------------------------------------------------------------ */
  for(xk = x, k = 0; k < n; xk += n, k++)
    d[k] = xk[k];
  for(j = 0; j < n; j++)
    perm[j] = j;
/* ------------------------------------------------------------
   Pivot in step k=0:n-1 on imax:
   ------------------------------------------------------------ */
  for(k = 0; k < n; k++){
/* ------------------------------------------------------------
   Let [imax,dk] = max(d(k:m))
   ------------------------------------------------------------ */
    dk = d[k]; imax = k;
    for(i = k + 1; i < n; i++)
      if(d[i] > dk){
        imax = i;
        dk = d[i];
      }
/* ------------------------------------------------------------
   k-th pivot is j=perm[imax].
   ------------------------------------------------------------ */
    d[imax] = d[k];
    j = perm[imax];                     /* original node number */
    uk = u + j * n;
    rowuk = u + k;
    perm[imax] = perm[k];
    perm[k] = j;
/* ------------------------------------------------------------
    Let uk[k] = (dk := sqrt(dk))
    ------------------------------------------------------------ */
    if(dk <= 0.0)
      return 1;      /* Matrix is not positive definite */
    dk = sqrt(dk);
    uk[k] = dk;
/* ------------------------------------------------------------
   Let  u(k,:) = (x(k,:)-uk'*u(0:k-1,:)) / dk, then d(k+1:n) -= u(k,:).^2.
   ------------------------------------------------------------ */
    for(i = k + 1; i < n; i++){
      icol = perm[i] * n;
      uki = (x[icol + j] - realdot(uk, u+icol, k)) / dk;
      rowuk[icol] = uki;
      d[i] -= SQR(uki);              /* d(:) -= u(k,:).^2 */
    }
  }
  return 0; /* success */
}

/* ************************************************************
   PROCEDURE prpicholpivot - U'*U factorization for nxn matrix,
     with column pivotting. Yields U with 1 >= \|U(i,i+1:n)\| forall i.
   INPUT
     n - order of matrix to be factored
     x,xpi - Full nxn. Matrix to be factored.
   OUTPUT
     u,upi - Full nxn, Cholesky factor: X=U'*U with U(:,perm)
       upper triangular.
     perm - Column permutation for maximal stability.
   WORK
     d - length n vector; the positive diagonal diag(U(:,perm)).^2,
         has decreasing order.
   RETURNS:
     0 = "success", 1 = "X is NOT positive definite"
   ************************************************************ */
int prpicholpivot(double *u, double *upi, int *perm, const double *x,
                  const double *xpi, const int n, double *d)
{
  int i,j,k,imax, icol;
  double *uk, *rowuk, *ukpi, *rowukpi;
  double dk, uki;
  const double *xk;
/* ------------------------------------------------------------
   Initialize: d = diag(x), perm = 0:n-1.
   ------------------------------------------------------------ */
  for(xk = x, k = 0; k < n; xk += n, k++)
    d[k] = xk[k];
  for(j = 0; j < n; j++)
    perm[j] = j;
/* ------------------------------------------------------------
   Pivot in step k=0:n-1 on imax:
   ------------------------------------------------------------ */
  for(k = 0; k < n; k++){
/* ------------------------------------------------------------
   Let [imax,dk] = max(d(k:m))
   ------------------------------------------------------------ */
    dk = d[k]; imax = k;
    for(i = k + 1; i < n; i++)
      if(d[i] > dk){
        imax = i;
        dk = d[i];
      }
/* ------------------------------------------------------------
   k-th pivot is j=perm[imax].
   ------------------------------------------------------------ */
    d[imax] = d[k];
    j = perm[imax];                     /* original node number */
    uk = u + j * n;
    rowuk = u + k;
    ukpi = upi + j * n;
    rowukpi = upi + k;
    perm[imax] = perm[k];
    perm[k] = j;
/* ------------------------------------------------------------
    Let uk[k] = (dk := sqrt(dk))
    ------------------------------------------------------------ */
    if(dk <= 0.0)
      return 1;      /* Matrix is not positive definite */
    dk = sqrt(dk);
    uk[k] = dk;
/* ------------------------------------------------------------
   Let  u(k,:) = (x(k,:)-uk'*u(0:k-1,:)) / dk,
   then d(k+1:n) -= u(k,:).^2.
   ------------------------------------------------------------ */
    for(i = k + 1; i < n; i++){
      icol = perm[i] * n;
      uki = x[icol + j];                                 /* real part */
      uki -= realdot(uk, u+icol, k) + realdot(ukpi, upi+icol,k);
      rowuk[icol] = (uki /= dk);
      d[i] -= SQR(uki);                               /* d(:) -= u(k,:).^2 */
      uki = xpi[icol + j];                               /* imaginary part */
      uki += realdot(ukpi, u+icol, k) - realdot(uk, upi+icol,k);
      rowukpi[icol] = (uki /= dk);
      d[i] -= SQR(uki);                               /* d(:) -= u(k,:).^2 */
    }
  }
  return 0; /* success */
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   [qdetx,ux,ispos,perm] = factorK(x,K);
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  mxArray *myplhs[NPAROUT];
  coneK cK;
  int i,k,nk,nksqr, sdplen,sdpdim,lenfull, fwsiz, ispos;
  const double *x;
  double *ux, *fwork, *permPr, *qdetx, *up, *uppi;
  int *iwork, *perm;
  double uxk;
  char use_pivot;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "factorK requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "factorK produces less output arguments");
  use_pivot = (nlhs == NPAROUT);
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
/* ------------------------------------------------------------
   Compute statistics: sdpdim = rdim+hdim, sdplen = sum(K.s).
   ------------------------------------------------------------ */
  lenfull = cK.lpN +  cK.qDim + cK.rDim + cK.hDim;
  sdpdim = cK.rDim + cK.hDim;
  sdplen = cK.rLen + cK.hLen;
/* ------------------------------------------------------------
   Get input vector x, skip LP part
   ------------------------------------------------------------ */
  mxAssert(mxGetM(X_IN) * mxGetN(X_IN) == lenfull, "x size mismatch.");
  x = mxGetPr(X_IN) + cK.lpN;
/* ------------------------------------------------------------
   Allocate output qdetx(lorN), UX(sdpdim), perm(sdplen), ispos(1).
   ------------------------------------------------------------ */
  QDETX_OUT = mxCreateDoubleMatrix(cK.lorN, 1, mxREAL);
  qdetx = mxGetPr(QDETX_OUT);
  UX_OUT = mxCreateDoubleMatrix(sdpdim, 1, mxREAL);
  ux = mxGetPr(UX_OUT);
  ISPOS_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
  PERM_OUT =  mxCreateDoubleMatrix(sdplen, 1, mxREAL);
  permPr = mxGetPr(PERM_OUT);
/* ------------------------------------------------------------
   Allocate working arrays iwork(sdplen),
   fwork(MAX(rmaxn^2,2*hmaxn^2) + MAX(rmaxn,hmaxn))
   ------------------------------------------------------------ */
  iwork = (int *) mxCalloc(sdplen, sizeof(int));
  perm = iwork;
  fwsiz = MAX(cK.rMaxn,cK.hMaxn);
  fwork = (double *) mxCalloc(fwsiz + MAX(SQR(cK.rMaxn),2*SQR(cK.hMaxn)),
                              sizeof(double));
  up = fwork + fwsiz;
  uppi = up + SQR(cK.hMaxn);
/* ------------------------------------------------------------
   LORENTZ:  qdetx = sqrt(qdet(x))
   ------------------------------------------------------------ */
  ispos = 1;
  for(k = 0; k < cK.lorN; k++){
    nk = cK.lorNL[k];
    if( (uxk = qdet(x,nk)) < 0.0){
      ispos = 0;
      break;
    }
    else
      qdetx[k] = sqrt(uxk);
    x += nk;
  }
/* ------------------------------------------------------------
   PSD: Cholesky factorization. If use_pivot, then do pivoting.
   ------------------------------------------------------------ */
  if(use_pivot){
    if(ispos)
      for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
        nk = cK.sdpNL[k];
        if(cholpivot(up,perm, x,nk, fwork)){
          ispos = 0;
          break;
        }
        uperm(ux, up, perm, nk);
        triu2sym(ux,nk);
        nksqr = SQR(nk);
        x += nksqr; ux += nksqr;
        perm += nk;
      }
/* ------------------------------------------------------------
   Complex Hermitian PSD pivoted Cholesky factorization
   ------------------------------------------------------------ */
    if(ispos)
      for(; k < cK.sdpN; k++){                    /* complex Hermitian */
        nk = cK.sdpNL[k];
        nksqr = SQR(nk);
        if(prpicholpivot(up,uppi,perm, x,x+nksqr,nk, fwork)){
          ispos = 0;
          break;
        }
        uperm(ux, up, perm, nk);                  /* real part */
        uperm(ux+nksqr, uppi, perm, nk);          /* imaginary part */
        triu2herm(ux,ux+nksqr,nk);
        nksqr += nksqr;                           /* 2*n^2 for real+imag */
        x += nksqr; ux += nksqr;
        perm += nk;
      }
/* ------------------------------------------------------------
   Convert "perm" to Fortran-index in doubles.
   ------------------------------------------------------------ */
    for(i = 0; i < sdplen; i++)
      permPr[i] = 1.0 + iwork[i];
  }
/* ------------------------------------------------------------
   PSD, !use_pivot: Cholesky without pivoting.
   First let ux = x, then ux=chol(ux).
   ------------------------------------------------------------ */
  else{           /* Cholesky real sym PSD without pivoting */
    if(ispos){
      memcpy(ux, x, sdpdim * sizeof(double));       /* copy real + complex */
      for(k = 0; k < cK.rsdpN; k++){                /* real symmetric */
        nk = cK.sdpNL[k];
        if(cholnopiv(ux,nk)){
          ispos = 0;
          break;
        }
        triu2sym(ux,nk);
        ux += SQR(nk);
      }
    }
/* ------------------------------------------------------------
   Complex Hermitian PSD Cholesky factorization, no pivoting.
   ------------------------------------------------------------ */
    if(ispos)
      for(; k < cK.sdpN; k++){                    /* complex Hermitian */
        nk = cK.sdpNL[k];
        nksqr = SQR(nk);
        if(prpicholnopiv(ux,ux+nksqr,nk)){
          ispos = 0;
         break;
        }
        triu2herm(ux,ux+nksqr,nk);
        ux += 2 * nksqr;
      }
  } /* !use_pivot */
/* ------------------------------------------------------------
   Return parameter ispos
   ------------------------------------------------------------ */
  *mxGetPr(ISPOS_OUT) = ispos;
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(iwork);
  mxFree(fwork);
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  i = MAX(nlhs, 1);
  memcpy(plhs,myplhs, i * sizeof(mxArray *));
  for(; i < NPAROUT; i++)
    mxDestroyArray(myplhs[i]);
}
