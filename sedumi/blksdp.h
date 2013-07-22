/* ************************************************************
   HEADER blksdp.h
   For use with mex-files in self-dual-minimization package.

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

#if !defined(BLKSDP)
#define BLKSDP
#include "mex.h"

/* ------------------------------------------------------------
   Type definitions:
   ------------------------------------------------------------ */
typedef struct{
 double *pr;
 mwIndex *jc, *ir;
    } jcir;

typedef struct{
 double *pr;
 const mwIndex *jc;
 const mwIndex *ir;
    } constjcir;

typedef struct{
  mwSize frN,lpN,lorN,rconeN,sdpN, rsdpN;
  mwSize qMaxn,rMaxn,hMaxn, rLen,hLen,  qDim,rDim,hDim;
  const double *lorNL,*rconeNL,*sdpNL;
} coneK;

/* ------------------------------------------------------------
   Macros:
   ------------------------------------------------------------ */
#if !defined(SQR)
#define SQR(x) ((x)*(x))
#endif

#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define  MIN(A, B)   ((A) < (B) ? (A) : (B))
#endif

#if !defined(SIGN)
#define  SIGN(A)   (2 * ((A) >= 0) - 1)
#endif

#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880	/* sqrt(2) */
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2   0.70710678118654752440	/* 1/sqrt(2) */
#endif

/* ************************************************************
   INT COMPARE: for searching an mwIndex array
   NOTE: qsort sorts in ascending (0,1,2,..) order, if the compare
     function returns  < 0 iff a<b, 0 iff a==b, > 0 iff a > b.
   ************************************************************ */
#if !defined(_COMPFUN_)
#define _COMPFUN_
typedef int (*COMPFUN)(const void *pa,const void *pb);
#endif

#define ibsearch(key,vec,n)  bsearch((void *)(key), (void *)(vec), (n), sizeof(mwSize), (COMPFUN) icmp)
#define iqsort(vec,n) qsort((void *)(vec), (n), sizeof(mwSize), (COMPFUN) icmp)

/* --------------------------------------
   KEY COMPARE: FOR SORTING AN (INT or FLOAT) ARRAY WITH INT-KEYS.
   -------------------------------------- */
typedef struct{
  mwIndex i,k;
} keyint;

#if !defined(_KEYDOUBLE_)
#define _KEYDOUBLE
typedef struct{
  double r;
  mwIndex k;
} keydouble;
#endif

#define kiqsort(vec,n)  qsort((void *)(vec), (n), sizeof(keyint), (COMPFUN) kicmp);
#define kdsortdec(vec,n)  qsort((void *)(vec), (n), sizeof(keydouble), (COMPFUN) kdcmpdec);

/*BLAS functions returning anything other than void need to be declared here as
 *Matlab does not include a header file, so the compiler will
 *assume they return an int.*/
#ifdef PC
#define FORT(x) x
#else
#define FORT(x) x ## _
#endif

double FORT(ddot)(const mwIndex*, const double*, const mwIndex*, const double*, const mwIndex*);
double FORT(dnrm2)(const mwIndex*, const double*, const mwIndex*);
mwIndex FORT(idamax)(const mwIndex *,const double *,const mwIndex *);
void FORT(dcopy)(const mwIndex*,const double*,const mwIndex*,double*,const mwIndex*);
void FORT(dscal)(const mwIndex*,const double*,double*,const mwIndex*);
void FORT(daxpy)(const mwIndex*,const double*,const double*,const mwIndex*,double*,const mwIndex*);

/* ------------------------------------------------------------
   Prototypes:
   ------------------------------------------------------------ */
char icmp(const mwIndex *a, const mwIndex *b);
bool intbsearch(mwIndex *pi, const mwIndex *x, const mwIndex n, const mwIndex key);
char intmbsearch(mwIndex *z, bool *found, const mwIndex *x, const mwIndex xnnz,
		const mwIndex *y, const mwIndex ynnz, mwIndex *iwork, const mwIndex iwsize);
char kicmp(const keyint *a, const keyint *b);
char kdcmpdec(const keydouble *a, const keydouble *b);
double realssqr(const double *x, const mwSize n);
double realdot(const double *x, const double *y, const mwSize n);
double selrealdot(const double *x, const double *y,
		  const mwIndex *sel, const mwSize nnz);
double realdotrow(const double *x, const double *y, const mwSize n);
void fromto(mwIndex *x, mwSize i, const mwSize n);
double triudotprod(const double *x, const double *y, const mwSize n);
double striudotprod(const double *x, const double *y, const mwSize n);
void tril2sym(double *r, const mwSize n);
void tril2herm(double *r, double *rpi, const mwSize n);
void triu2sym(double *r, const mwSize n);
void triu2herm(double *r, double *rpi, const mwSize n);
void scalarmul(double *r, const double alpha,const double *x,const mwSize n);
void addscalarmul(double *r, const double alpha,const double *x,const mwSize n);
void subscalarmul(double *x, const double alpha, const double *y, const mwSize n);
void realHadamard(double * r, const double *x, const double *y, const mwSize n);
void minusHadamard(double * r, const double *x, const double *y, const mwSize n);
void realHadarow(double * r, const double *x, const double *y, const mwSize n);
void realHadadiv(double * r, const double *x, const double *y, const mwSize n);
void fzeros(double *z,const mwSize n);
void conepars(const mxArray *mxK, coneK *pK);
void someStats(mwSize *pxmax, mwIndex *pxsum, mwIndex *pxssqr,
	       const double *x, const mwSize n);
mwIndex spsqrscale(double *z, mwIndex *blks, const mwIndex *zjc, const mwIndex *zir,
               const mwIndex *znnz, const double *d,
               const mwIndex *xir, const double *xpr, mwIndex xjc0, const mwIndex xjc1,
               const mwIndex *blkstart, const mwIndex *xblk, const mwIndex *psdNL,
               const mwIndex rpsdN, double *fwork, mwIndex *iwork);
#ifdef OLDSEDUMI
double qscale(double *z,const double *x,const double *y,
              const double rdetx,const mwIndex n);
void qlmul(double *z,const double *x,const double *y,
	   const double rdetx,const mwIndex n);
void qldiv(double *z,const double *x,const double *y,
	   const double rdetx,const mwIndex n);
void vec2blks(mwIndex *blklocs, const mwIndex *blkstart, const mwIndex *yir,
              const mwIndex ystart, const mwIndex ynnz, const mwIndex nblk);
void vec2selblks(mwIndex *blklocs, const mwIndex *blkstart, const mwIndex *yir,
                 const mwIndex ystart, const mwIndex ynnz,
                 const mwIndex *blkir, const mwIndex blknnz);
mwIndex lqdsqrx(double *z,
            const mwIndex *xir, const double *xpr, const mwIndex xjc0,
            const mwIndex xjcq, const mwIndex xjcs, const mwIndex *qir,
            const mwIndex *blkstart,
            const double *dsqr, const double *detd);
mwIndex blkpsdscale(double *z, const mwIndex *zir, const mwIndex zjc1,
		const double *u, const mwIndex *invperm, const double *x,
		const mwIndex *xblk, const mwIndex blkjc0, const mwIndex blkjc1,
		const mwIndex *blkstart, const mwIndex *psdNL, const mwIndex *cumpsdNL,
		const mwIndex rpsdN, double *fwork);
#endif
void uperm(double *y, const double *u, const mwIndex *perm, const mwIndex n);
/* ------------------------------------------------------------
   For auxfwdpr1:
   ------------------------------------------------------------ */
void fwipr1(double *y, const double *p, const double *beta,
            const mwSize m, const mwSize n);
void fwipr1o(double *y, const mwIndex *perm, const double *p, const double *beta,
             const mwSize m, const mwSize n);
#endif
