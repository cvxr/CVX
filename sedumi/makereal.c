/* ************************************************************
%                                                y = makereal(x,K,cpx)
% MAKEREAL  Converts matrix in MATLAB-complex format to internal
%  SeDuMi format.
%
% SEE ALSO sedumi.
% ******************** INTERNAL FUNCTION OF SEDUMI ********************

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

#include <math.h>
#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define Y_OUT plhs[0]
#define NPAROUT 1

#define X_IN prhs[0]
#define K_IN prhs[1]
#define CPX_IN prhs[2]
#define NPARIN 3

/* ************************************************************
   PROCEDURE intscalaradd - Let y = d + x.
   ************************************************************ */
void intscalaradd(mwIndex *y, const mwSize d, const mwIndex *x, const mwSize n)
{
  mwIndex i;
  for(i = 0; i < n; i++)
    y[i] = x[i] + d;
}

/* ************************************************************
   PROCEDURE writenz
   ************************************************************ */
void writenz(mwIndex *yir, double *ypr, mwIndex *pjnz, const mwIndex *xir,
             const double *xpr, const mwIndex jc0, const mwIndex jc1, const mwIndex d)
{
  mwIndex i,jnz;
  double yi;
  jnz = *pjnz;
/* ------------------------------------------------------------
   Write nonzeros from x into y that are REALLY nonzero; also
   translate subscript by d.
   ------------------------------------------------------------ */
  for(i = jc0; i < jc1; i++)
    if( (yi = xpr[i]) != 0.0){          /* Really nonzero? */
      ypr[jnz] = yi;
      yir[jnz] = xir[i] + d;
      jnz++;
    }
  *pjnz = jnz;
}

/* ************************************************************
   PROCEDURE fmakereal
   INPUT
     xir,xpi - sparse imaginary vector with xnnz "nonzeros" in LP/Lorentz.
        If xpi == NULL, then all imaginary data is zero.
     xnnz - length of xir,xpi, before PSD part.
     cpxf - length nf integer array, listing free imaginary vars.
     nf - length of cpxf
     iwsize Length of iwork, should be 2+nf + floor(log_2(1+nf)).
   OUTPUT
     yir, ypr - sparse real output vector, ynnz nonzeros.
   WORK ARRAYS
     cfound - length MAXN := MAX(nf,nx,2*ns) character working array.
     iwork - lengt iwsize working array. Needs
      iwsize >= 2+nf + floor(log_2(1+nf)).
   RETURNS ynnz (<=2*xnnz), number of free imag nonzeros in y.
   ************************************************************ */
mwIndex fmakereal(mwIndex *yir, double *ypr, const mwIndex *xir, const double *xpi,
	     const mwIndex xnnz, const mwIndex *cpxf, const mwIndex nf,
	     bool *cfound, mwIndex *iwork, mwIndex iwsize)
{
  mwIndex i,jnz;
  mwIndex *ipos;
  double yj;

  if(xpi == (double *) NULL)
    return 0;                  /* No imaginary nonzeros */
/* ------------------------------------------------------------
   Partition WORKING ARRAY:
   ------------------------------------------------------------ */
  ipos = iwork;
  iwork += 2 + nf;
  iwsize -= 2 + nf;
/* ------------------------------------------------------------
   Locate cpx.f in xir(1:xnnz)
   ------------------------------------------------------------ */
  if(intmbsearch(ipos, cfound, xir, xnnz, cpxf, nf, iwork, iwsize) != 0)
    mexErrMsgTxt("Out of working space");
/* ------------------------------------------------------------
   Write y = sparse(imag(x(cpx.f))), the free imaginary components
   ------------------------------------------------------------ */
  jnz = 0;
  for(i = 0; i < nf; i++)
    if(cfound[i] != 0){
      if( (yj = xpi[ipos[i+1]]) != 0.0){
	ypr[jnz] = yj;
	yir[jnz] = i;
	jnz++;
      }
    }
/* ------------------------------------------------------------
   RETURN number of (free imag) nonzeros in y
   ------------------------------------------------------------ */
  return jnz;
}


/* ************************************************************
   PROCEDURE xmakereal
   INPUT
     xir,xpr,xpi - sparse vector with xnnz nonzeros in LP/Lor.
     xnnz - length of xir,xpr,xpi; counting only LP/Lor nonzeros.
     idelta - Subscript adjustment for x(1), i.e. nf.
     cpxx - length nx integer array, listing Lorentz constrained
       imaginary vars.
     nx - length of cpxx.
     iwsize Length of iwork, should be 2+nx + floor(log_2(1+nx)).
   OUTPUT
     yir, ypr - sparse real output vector, ynnz nonzeros.
   WORK ARRAYS
     cfound - length nx character working array.
     iwork - lengt iwsize working array. Needs
      iwsize >= 2+nx + floor(log_2(1+nx)).
   RETURNS ynnz (<=2*xnnz), number of nonzeros in y.
   ************************************************************ */
mwIndex xmakereal(mwIndex *yir, double *ypr,
	      const mwIndex *xir, const double *xpr, const double *xpi,
	      const mwIndex xnnz, const mwIndex idelta, const mwIndex *cpxx, const mwIndex nx,
	      bool *cfound, mwIndex *iwork, mwIndex iwsize)
{
  mwIndex i,j,k,inz,jnz;
  mwIndex *ipos;
  double yj;
/* ------------------------------------------------------------
   Partition WORKING ARRAY:
   ------------------------------------------------------------ */
  ipos = iwork;
  iwork += 2 + nx;
  iwsize -= 2 + nx;
/* ------------------------------------------------------------
   Locate cpx.x in xir(1:xnnz)
   ------------------------------------------------------------ */
  if(intmbsearch(ipos, cfound, xir, xnnz, cpxx, nx, iwork, iwsize) != 0)
    mexErrMsgTxt("Out of working space");
/* ------------------------------------------------------------
   Write y(nf:nf+dimflqr+nx) = merge( real(x), imag(x(cpx.xcomplex)) )
   ------------------------------------------------------------ */
  memcpy(ypr, xpr, ipos[1] * sizeof(double));
  if(idelta != 0)
    intscalaradd(yir, idelta, xir, ipos[1]);
  else
    memcpy(yir, xir, ipos[1] * sizeof(mwIndex));
  jnz = ipos[1];
  for(i = 0; i < nx; ){
    if(cfound[i++]){
      inz = ipos[i];
      j = xir[inz] + idelta + i;         /* new index of imag part */
      if((yj = xpr[inz]) != 0.0){
	ypr[jnz] = yj;
	yir[jnz] = j-1;              /* index of real part */
	jnz++;
      }
      if(xpi != (double *) NULL)
        if((yj = xpi[inz]) != 0.0){
          ypr[jnz] = yj;
          yir[jnz] = j;                /* imag part */
          jnz++;
        }
      inz++;   /* point to first nonzero in x beyond cpx.x(i) */
    }
    else
      inz = ipos[i];    /*  point to first nonzero in x beyond cpx.x(i) */
/* ------------------------------------------------------------
   Let y(nf+i+(cpx.x(i)+1:cpx.x(i+1)-1)) = real(x(cpx.x(i)+1:cpx.x(i+1)-1))
   If i==nx, then this is simply the remainder.
   ------------------------------------------------------------ */
    k = ipos[i+1] - inz;                  /* number of nonzeros to copy */
    memcpy(ypr + jnz, xpr + inz, k * sizeof(double));
    intscalaradd(yir+jnz, idelta + i, xir+inz, k);     /* adjust subscript */
    jnz += k;
  }
/* ------------------------------------------------------------
   RETURN number of nonzeros written into y
   ------------------------------------------------------------ */
  return jnz;
}

/* ************************************************************
   PROCEDURE smakereal
   INPUT
     xir,xpr,xpi - sparse vector with xnnz nonzeros, in PSD ONLY!.
     xnnz - length of xir,xpr,xpi (PSD only).
     idelta - subscript adjustment of 1st PSD entry
     cpxsi - length twons increasing integer array. The old subscripts
       of the kth Hermitian block in x are cpxsi[2*k]:cpxsi[2*k+1]-1.
     twons - twons/2 is number of Hermitian PSD blocks.
     lenfull - length of full(x)-vector, 1 beyond last possible subscript.
     iwsize Length of iwork, should be 2*(1+ns) + floor(log_2(1+2*ns)).
   OUTPUT
     yir, ypr - sparse real output vector, ynnz nonzeros.
   WORK ARRAYS
     cfound - length MAXN := MAX(nf,nx,2*ns) character working array.
     iwork - lengt iwsize working array. Needs
      iwsize >= 2+2*ns + floor(log_2(1+2*ns)).
   RETURNS ynnz (<=2*xnnz), number of nonzeros in y.
   ************************************************************ */
mwIndex smakereal(mwIndex *yir, double *ypr,
	      const mwIndex *xir, const double *xpr, const double *xpi,
	      const mwIndex xnnz, const mwIndex idelta,
	      const mwIndex *cpxsi, const mwIndex twons, const mwIndex lenfull,
	      bool *cfound, mwIndex *iwork, mwIndex iwsize)
{
  mwIndex i,j,k,inz,jnz;
  mwIndex *ipos;
/* ------------------------------------------------------------
   Partition WORKING ARRAY:
   ------------------------------------------------------------ */
  ipos = iwork;         /* length 2*(1+ns) array */
  iwork += 2 + twons;
  iwsize -= 2 + twons;
/* ------------------------------------------------------------
   Locate cpxsi in x(jcs:end); these mark the start+ends of Hermitian blocks,
   and hence the complement are real blocks.
   ------------------------------------------------------------ */
  if(intmbsearch(ipos, cfound, xir, xnnz, cpxsi, twons, iwork, iwsize) != 0)
    mexErrMsgTxt("Out of working space");
/* ------------------------------------------------------------
   Write real PSD blocks into y, i.e. skip Hermitian blocks
   ------------------------------------------------------------ */
  memcpy(ypr, xpr, ipos[1] * sizeof(double));
  if(idelta != 0)
    intscalaradd(yir, idelta, xir, ipos[1]);
  else
    memcpy(yir, xir, ipos[1] * sizeof(mwIndex));
  jnz = ipos[1];
  j = idelta;             /* subscript adjustment */
  for(i = 0; i < twons; ){
    j -= cpxsi[i+1] - cpxsi[i];               /* skip complex block */
    i += 2;
    inz = ipos[i];
    k = ipos[i + 1] - inz;
    memcpy(ypr + jnz, xpr + inz, k * sizeof(double));
    intscalaradd(yir+jnz, j, xir+inz, k);
    jnz += k;
  }
/* ------------------------------------------------------------
   Write Hermitian PSD blocks into y
   ------------------------------------------------------------ */
  j += lenfull;    /* j points to 1st available index for Hermitian blocks */
  for(i = 0; i < twons; i += 2){
    k = cpxsi[i];              /* Old 1st index of Herm PSD block */
/* ---------- write real part ---------- */
    writenz(yir, ypr, &jnz, xir, xpr, ipos[i+1],ipos[i+2],j-k);
    j += cpxsi[i+1] - k;                 /* point to 1st available index */
/* ---------- write imag part ---------- */
    if(xpi != (double *) NULL)
      writenz(yir, ypr, &jnz, xir, xpi, ipos[i+1],ipos[i+2],j-k);
    j += cpxsi[i+1] - k;                 /* point to 1st available index */
  }
/* ------------------------------------------------------------
   RETURN number of nonzeros written into y
   ------------------------------------------------------------ */
  return jnz;
}

/* ************************************************************
   PROCEDURE makereal
   INPUT
     xir,xpr,xpi - sparse vector with xnnz nonzeros.
     xnnz - length of xir,xpr,xpi.
     cpxf - length nf integer array, listing free imaginary vars.
     nf - length of cpxf
     cpxx - length nx integer array, listing Lorentz constrained
       imaginary vars.
     nx - length of cpxx.
     cpxsi - length 2*ns increasing integer array. The old subscripts
       of the kth Hermitian block in x are cpxsi[2*k]:cpxsi[2*k+1]-1.
     ns - number of Hermitian PSD blocks.
     lenfull - length of full(x)-vector, 1 beyond last possible subscript.
     iwsize Length of iwork, should be 2 + MAXN + floor(log_2(1+MAXN)).
     lenfull - dimension of x.
     dimflqr - dimension of f/l/q/r part in x, i.e. 1st valid PSD subscript.
   OUTPUT
     yir, ypr - sparse real output vector, ynnz nonzeros.
   WORK ARRAYS
     cfound - length MAXN := MAX(nf,nx,2*ns) character working array.
     iwork - lengt iwsize working array. Needs
      iwsize >= 2 + MAXN + floor(log_2(1+MAXN)).
   RETURNS ynnz (<=2*xnnz), number of nonzeros in y.
   ************************************************************ */
mwIndex makereal(mwIndex *yir, double *ypr,
             const mwIndex *xir, const double *xpr, const double *xpi,
             const mwIndex xnnz, const mwIndex *cpxf, const mwIndex nf,
             const mwIndex *cpxx, const mwIndex nx, const mwIndex *cpxsi,
             const mwIndex ns, const mwIndex lenfull, const mwIndex dimflqr,
             bool *cfound, mwIndex *iwork, mwIndex iwsize)
{
  mwIndex jcs, jnz;
/* ------------------------------------------------------------
   Find position of 1st PSD nonzero
   ------------------------------------------------------------ */
  jcs = 0;
  intbsearch(&jcs, xir, xnnz, dimflqr);
/* ------------------------------------------------------------
   Write free imaginary nonzeros into y
   ------------------------------------------------------------ */
  jnz = fmakereal(yir,ypr, xir,xpi,jcs, cpxf,nf, cfound,iwork,iwsize);
/* ------------------------------------------------------------
   Write LP/Lorentz nonzeros into y, make Lorentz bounded imag nonzeros
   explicit reals.
   ------------------------------------------------------------ */
  jnz += xmakereal(yir+jnz,ypr+jnz, xir,xpr,xpi,jcs, nf, cpxx,nx,
                   cfound,iwork,iwsize);
/* ------------------------------------------------------------
   Write PSD nonzeros into y. First the real blocks, then Hermitian.
   ------------------------------------------------------------ */
  if(xpi != (double *) NULL)
    xpi += jcs;                   /* point to 1st imag PSD nonzero */
  jnz += smakereal(yir+jnz, ypr+jnz, xir+jcs,xpr+jcs,xpi,xnnz-jcs,
                   nf+nx, cpxsi,2*ns, lenfull, cfound,iwork,iwsize);
/* ------------------------------------------------------------
   Return number of nonzeros in y
   ------------------------------------------------------------ */
  return jnz;
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
  mwIndex i,j,jnz,MAXN, m,dimflqr,lenfull,reallength, nf,nx,ns, iwsize;
  mwIndex *iwork, *sdpNL, *cpxf, *cpxx, *cpxs, *cpxsi;
  bool *cwork;
  const double *cpxfPr, *cpxxPr, *cpxsPr, *xpi;
  mxArray *MY_FIELD;
  coneK cK;
  jcir x,y;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "makereal requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "makereal produces less output arguments");
/* ------------------------------------------------------------
   Disassemble cone K structure
   ------------------------------------------------------------ */
  conepars(K_IN, &cK);
  dimflqr = cK.frN + cK.lpN + cK.qDim;
  for(i = 0; i < cK.rconeN; i++)        /* add dim of rotated cone */
    dimflqr += cK.rconeNL[i];
  lenfull =  dimflqr + cK.rDim + cK.hDim;
/* ------------------------------------------------------------
   Disassemble cpx structure
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(CPX_IN), "Parameter `cpx' should be a structure.");
  if( (MY_FIELD = mxGetField(CPX_IN,(mwIndex)0,"f")) == NULL)  /* cpx.f */
    nf = 0;
  else{
    nf = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
    cpxfPr = mxGetPr(MY_FIELD);
  }
  if( (MY_FIELD = mxGetField(CPX_IN,(mwIndex)0,"s")) == NULL)  /* cpx.s */
    ns = 0;
  else{
    ns = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
    cpxsPr = mxGetPr(MY_FIELD);
  }
  if( (MY_FIELD = mxGetField(CPX_IN,(mwIndex)0,"x")) == NULL)  /* cpx.x */
    nx = 0;
  else{
    nx = mxGetM(MY_FIELD) * mxGetN(MY_FIELD);
    cpxxPr = mxGetPr(MY_FIELD);
  }
/* ------------------------------------------------------------
   Get input matrix x
   ------------------------------------------------------------ */
  mxAssert(mxGetM(X_IN) == lenfull, "Size x mismatch.");
  m = mxGetN(X_IN);                       /* number of columns to handle */
  mxAssert( mxIsSparse(X_IN), "X should be sparse.");
  x.pr = mxGetPr(X_IN);
  if(mxIsComplex(X_IN))
    xpi = mxGetPi(X_IN);
  else
    xpi = (double *) NULL;
  x.jc = mxGetJc(X_IN);
  x.ir = mxGetIr(X_IN);
/* ------------------------------------------------------------
   Allocate iwork[iwsiz], cwork[MAXN], {K.s, cpx.{f[nf],x[nx],s[ns],si[2*ns]}}
   ------------------------------------------------------------ */
  MAXN = MAX(MAX(nf,nx),2*ns);
  iwsize = floor(log(1.0 + MAXN) / log(2.0));       /* for binary tree search */
  iwsize += 2 * ns + 2 + MAXN;
  iwork = (mwIndex *) mxCalloc(iwsize, sizeof(mwIndex));
  cwork = (bool *) mxCalloc(MAX(1,MAXN), sizeof(bool));
  sdpNL = (mwIndex *) mxCalloc(MAX(1, cK.sdpN + nf + nx + 3*ns), sizeof(mwIndex));
  cpxf = sdpNL + cK.sdpN;
  cpxx = cpxf + nf;
  cpxs = cpxx + nx;
  cpxsi = cpxs + ns;        /* length 2*ns */
/* ------------------------------------------------------------
   Convert double to mwIndex
   ------------------------------------------------------------ */
  for(i = 0; i < cK.sdpN; i++){        /* K.s */
    j = cK.sdpNL[i];
    sdpNL[i] = j;                    /* These are lengths, not subscripts!*/
  }
  for(i = 0; i < nf; i++){             /* cpx.f */
    j = cpxfPr[i];
    cpxf[i] = --j;
  }
  for(i = 0; i < nx; i++){             /* cpx.x */
    j = cpxxPr[i];
    cpxx[i] = --j;
  }
  for(i = 0; i < ns; i++){             /* cpx.s */
    j = cpxsPr[i];
    cpxs[i] = --j;
  }
/* ------------------------------------------------------------
   Create cpxsi(1:2*ns). This lists the 1st subscript and the
   1-beyond-last subscript of each Hermitian PSD matrix in x.
   ------------------------------------------------------------ */
  jnz = dimflqr;
  j = 0;
  for(i = 0; i < ns; i++){
    for(; j < cpxs[i]; j++)
      jnz += SQR(sdpNL[j]);
    cpxsi[2 * i] = jnz;            /* start of Hermitian block */
    jnz += SQR(sdpNL[j]);
    j++;
    cpxsi[2 * i + 1] = jnz;        /* end of Hermitian block */
  }
/* ------------------------------------------------------------
   Allocate output Y = sparse([],[],[],reallength(x),m,2 * nnz(x))
   ------------------------------------------------------------ */
  reallength = lenfull + nf + nx;
  for(i = 0; i < ns; i++){
    reallength += cpxsi[2*i+1] - cpxsi[2*i];
  }
  jnz = x.jc[m];
  if(xpi != (double *) NULL)
    jnz *= 2;                       /* reserve room for imaginary parts */
  Y_OUT = mxCreateSparse(reallength, m, jnz, mxREAL);
  y.pr = mxGetPr(Y_OUT);
  y.jc = mxGetJc(Y_OUT);
  y.ir = mxGetIr(Y_OUT);
/* ------------------------------------------------------------
   The real job, MAKEREAL:
   ------------------------------------------------------------ */
  y.jc[0] = 0;
  jnz = 0;
  for(i = 0; i < m; i++){
    y.jc[i] = jnz;
    j = x.jc[i+1] - x.jc[i];
    jnz += makereal(y.ir+jnz, y.pr+jnz, x.ir+x.jc[i],x.pr+x.jc[i],xpi,j,
                    cpxf,nf, cpxx,nx, cpxsi,ns, lenfull, dimflqr,
                    cwork, iwork, iwsize);
    if(xpi != (double *) NULL)
      xpi += j;                 /* To next imaginary column */
  }
  y.jc[i] = jnz;
  mxAssert(jnz <= mxGetNzmax(Y_OUT),"");
/* ------------------------------------------------------------
   REALLOC: Shrink Y to its current size
   ------------------------------------------------------------ */
  jnz = MAX(jnz,1);
  if( (y.pr = (double *) mxRealloc(y.pr, jnz*sizeof(double))) == NULL)
    mexErrMsgTxt("Memory reallocation error");
  mxSetPr(Y_OUT,y.pr);
  if( (y.ir = (mwIndex *) mxRealloc(y.ir, jnz*sizeof(mwIndex))) == NULL)
    mexErrMsgTxt("Memory reallocation error");
  mxSetIr(Y_OUT,y.ir);
  mxSetNzmax(Y_OUT,jnz);
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(sdpNL);
  mxFree(cwork);
  mxFree(iwork);
}
