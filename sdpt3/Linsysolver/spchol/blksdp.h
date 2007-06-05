/* ************************************************************
   HEADER blksdp.h
   For use with mex-files in self-dual-minimization package.

 %  
 %   This file is part of SeDuMi 1.02   (03AUG1998)
 %   Copyright (C) 1998 Jos F. Sturm
 %   CRL, McMaster University, Canada.
 %   Supported by the Netherlands Organization for Scientific Research (NWO).
 % 
 %   This program is free software; you can redistribute it and/or modify
 %   it under the terms of the GNU General Public License as published by
 %   the Free Software Foundation; either version 2 of the License, or
 %   (at your option) any later version.
 % 
 %   This program is distributed in the hope that it will be useful,
 %   but WITHOUT ANY WARRANTY; without even the implied warranty of
 %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 %   GNU General Public License for more details.
 % 
 %   You should have received a copy of the GNU General Public License
 %   along with this program; if not, write to the Free Software
 %   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 %

   ************************************************************ */
#if !defined(BLKSDP)
#define BLKSDP
#include "mex.h"

/* ------------------------------------------------------------
   Type definitions:
   ------------------------------------------------------------ */
typedef struct{
 double *pr;
 int *jc, *ir;
    } jcir;

typedef struct{
 double *pr;
 const int *jc;
 const int *ir;
    } constjcir;

typedef struct{
  int frN,lpN,lorN,rconeN,sdpN, rsdpN;
  int qMaxn,rMaxn,hMaxn, rLen,hLen,  qDim,rDim,hDim;
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
   INT COMPARE: for searching an int array
   NOTE: qsort sorts in ascending (0,1,2,..) order, if the compare
     function returns  < 0 iff a<b, 0 iff a==b, > 0 iff a > b.
   ************************************************************ */
#if !defined(_COMPFUN_)
#define _COMPFUN_
typedef int (*COMPFUN)(const void *pa,const void *pb);
#endif

#define ibsearch(key,vec,n)  bsearch((void *)(key), (void *)(vec), (n), sizeof(int), (COMPFUN) icmp)
#define iqsort(vec,n) qsort((void *)(vec), (n), sizeof(int), (COMPFUN) icmp)

/* --------------------------------------
   KEY COMPARE: FOR SORTING AN (INT or FLOAT) ARRAY WITH INT-KEYS.
   -------------------------------------- */
typedef struct{
  int i,k;
} keyint;

#if !defined(_KEYDOUBLE_)
#define _KEYDOUBLE
typedef struct{
  double r;
  int k;
} keydouble;
#endif

#define kiqsort(vec,n)  qsort((void *)(vec), (n), sizeof(keyint), (COMPFUN) kicmp);
#define kdsortdec(vec,n)  qsort((void *)(vec), (n), sizeof(keydouble), (COMPFUN) kdcmpdec);


/* ------------------------------------------------------------
   Prototypes:
   ------------------------------------------------------------ */
int icmp(const int *a, const int *b);
int intbsearch(int *pi, const int *x, const int n, const int key);
int intmbsearch(int *z, char *found, const int *x, const int xnnz,
		const int *y, const int ynnz, int *iwork, const int iwsize);
int kicmp(const keyint *a, const keyint *b);
int kdcmpdec(const keydouble *a, const keydouble *b);
double realssqr(const double *x, const int n);
double realdot(const double *x, const double *y, const int n);
double selrealdot(const double *x, const double *y,
		  const int *sel, const int nnz);
double realdotrow(const double *x, const double *y, const int n);
void fromto(int *x, int i, const int n);
double triudotprod(const double *x, const double *y, const int n);
double striudotprod(const double *x, const double *y, const int n);
void tril2sym(double *r, const int n);
void tril2herm(double *r, double *rpi, const int n);
void triu2sym(double *r, const int n);
void triu2herm(double *r, double *rpi, const int n);
void scalarmul(double *r, const double alpha,const double *x,const int n);
void addscalarmul(double *r, const double alpha,const double *x,const int n);
void subscalarmul(double *x, const double alpha, const double *y, const int n);
void realHadamard(double * r, const double *x, const double *y, const int n);
void minusHadamard(double * r, const double *x, const double *y, const int n);
void realHadarow(double * r, const double *x, const double *y, const int n);
void realHadadiv(double * r, const double *x, const double *y, const int n);
void fzeros(double *z,const int n);
void conepars(const mxArray *mxK, coneK *pK);
void someStats(int *pxmax, int *pxsum, int *pxssqr,
	       const double *x, const int n);
int spsqrscale(double *z, int *blks, const int *zjc, const int *zir,
               const int *znnz, const double *d,
               const int *xir, const double *xpr, int xjc0, const int xjc1,
               const int *blkstart, const int *xblk, const int *psdNL,
               const int rpsdN, double *fwork, int *iwork);
#ifdef OLDSEDUMI
double qscale(double *z,const double *x,const double *y,
              const double rdetx,const int n);
void qlmul(double *z,const double *x,const double *y,
	   const double rdetx,const int n);
void qldiv(double *z,const double *x,const double *y,
	   const double rdetx,const int n);
void vec2blks(int *blklocs, const int *blkstart, const int *yir,
              const int ystart, const int ynnz, const int nblk);
void vec2selblks(int *blklocs, const int *blkstart, const int *yir,
                 const int ystart, const int ynnz,
                 const int *blkir, const int blknnz);
int lqdsqrx(double *z,
            const int *xir, const double *xpr, const int xjc0,
            const int xjcq, const int xjcs, const int *qir,
            const int *blkstart,
            const double *dsqr, const double *detd);
int blkpsdscale(double *z, const int *zir, const int zjc1,
		const double *u, const int *invperm, const double *x,
		const int *xblk, const int blkjc0, const int blkjc1,
		const int *blkstart, const int *psdNL, const int *cumpsdNL,
		const int rpsdN, double *fwork);
#endif
void uperm(double *y, const double *u, const int *perm, const int n);
/* ------------------------------------------------------------
   For auxfwdpr1:
   ------------------------------------------------------------ */
void fwipr1(double *y, const double *p, const double *beta,
            const int m, const int n);
void fwipr1o(double *y, const int *perm, const double *p, const double *beta,
             const int m, const int n);
#endif
