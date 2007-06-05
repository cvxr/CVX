/*
 L = symfctmex(X, perm)
   Computes sparse symbolic factor L.L, updated permutation L.perm,
   super-node partition L.xsuper.

   Invokes SPARSPAK-A (ANSI FORTRAN) RELEASE III,
   by Joseph Liu (UNIVERSITY OF WATERLOO).
*/
/*
    This file is part of SeDuMi 1.03BETA
    Copyright (C) 1999 Jos F. Sturm
    Dept. Quantitative Economics, Maastricht University, the Netherlands.
    Affiliations up to SeDuMi 1.02 (AUG1998):
      CRL, McMaster University, Canada.
      Supported by the Netherlands Organization for Scientific Research (NWO).
  
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
  
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
  
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "mex.h"

#define L_OUT plhs[0]
#define NPAROUT 1

#define X_IN prhs[0]
#define PERM_IN prhs[1]
#define NPARIN 2
#ifdef DO_BFINIT
#define CACHSZ_IN prhs[2]
#endif

#if !defined(SQR)
#define SQR(x) ((x)*(x))
#endif

#if !defined(MIN)
#define  MIN(A, B)   ((A) < (B) ? (A) : (B))
#endif
#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif


/* ============================================================
   SUBROUTINES:
   ============================================================ */
/* ------------------------------------------------------------
   GETADJ - Copies off-diagonal entries from C-style sparse
     matrix (cjc,cir) to Fortran style sparse matrix (forjc,forir).
     On input, n is number of columns.
   ------------------------------------------------------------ */
void getadj(int *forjc,int *forir,const int *cjc,const int *cir,
            const int n)
{
  int i,j,inz,ix;
  
  inz = 0;
  for(j = 0; j < n; j++){
    forjc[j] = inz + 1;
    for(ix = cjc[j]; ix < cjc[j+1]; ix++)
      if((i = cir[ix]) != j)
        forir[inz++] = ++i;
  }
  forjc[n] = ++inz;
}

/* ------------------------------------------------------------
   EXPANDSUB -
   ------------------------------------------------------------ */
void expandsub(const int n,const int nsuper,
               const int *xsuper,const int *xlindx,
               int *Ljc,int *Lir)
{
  int j, jsup, jcol, ix, jpnt,ipnt;
/* ------------------------------------------------------------
   Convert Ljc from FORTRAN to C, i.e. -=1
   ------------------------------------------------------------ */
  for(j = 0; j <= n; j++)
    Ljc[j]-=1;
/* ------------------------------------------------------------
   For each snode: bring subscript to first column of snode,
   and translate from Fortran to C, i.e. -=1.
   ------------------------------------------------------------ */
  for(jsup = nsuper; jsup > 0; jsup--){
    jcol = xsuper[jsup-1];
    jpnt = Ljc[jcol];          /* points behind 1st column */
    ipnt = jpnt;
    for(ix = xlindx[jsup] - 1; ix >= xlindx[jsup-1]; )
      Lir[--ipnt] = Lir[--ix] - 1;
    if(ipnt != Ljc[jcol-1])
      mexErrMsgTxt("Input error expandsub.");
/* ------------------------------------------------------------
   Fill in subscripts of other columns in snode
   ------------------------------------------------------------ */
    for(; jcol < xsuper[jsup] - 1; jcol++){
      ipnt = jpnt;           /* behind 1st column */
      for(ix = Ljc[jcol+1]; ix > Ljc[jcol];)
        Lir[--ix] = Lir[--ipnt];
    }
  }
}

#define NL_FIELDS 3
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   L = symfctmex(X, perm, cachsz)
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  int m, i,j, iwsiz, nsuper, nsub, flag, nnzl;
  double *permPr,*xsuperPr,*Lpr, *outflag;
  int *Ljc, *Lir, *xadj, *adjncy, *Xjc, *Xir,
    *perm, *snode, *xsuper, *iwork,*xlindx,
    *invp, *colcnt;
  mxArray *L_FIELD;
  const char *LFieldnames[] = {"L", "perm", "xsuper"};
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARIN)
    mexErrMsgTxt("symfctmex requires more input arguments");
  if(nlhs > 2)
    mexErrMsgTxt("symfctmex produces less output arguments");
/* ------------------------------------------------------------
   Check input sizes (ADJ,perm)
   ------------------------------------------------------------ */
  if( (m = mxGetM(X_IN)) != mxGetN(X_IN) )
    mexErrMsgTxt("X must be square");
  if(!mxIsSparse(X_IN))
    mexErrMsgTxt("X must be sparse");
  if( (mxGetM(PERM_IN) * mxGetN(PERM_IN)) != m)
    mexErrMsgTxt("perm size mismatch");
/* ------------------------------------------------------------
   Get input (X,perm)
   ------------------------------------------------------------ */
  Xjc = mxGetJc(X_IN);
  Xir = mxGetIr(X_IN);
  permPr = mxGetPr(PERM_IN);
#ifdef DO_BFINIT
  cachsz = mxGetScalar(CACHSZ_IN);
#endif
/* ------------------------------------------------------------
   Allocate working arrays:
    int xadj(m+1), adjncy(Xnnz), perm(m), invp(m), colcnt(m), snode(m),
       xsuper(m+1), iwork(iwsize), xlindx(m+1), split(m)
   ------------------------------------------------------------ */
  xadj = (int *) mxCalloc(m+1,sizeof(int));
  adjncy = (int *) mxCalloc(Xjc[m], sizeof(int));
  perm   = (int *) mxCalloc(m,  sizeof(int));
  invp   = (int *) mxCalloc(m,  sizeof(int));
  colcnt = (int *) mxCalloc(m,  sizeof(int));
  snode  = (int *) mxCalloc(m,  sizeof(int));
  xsuper = (int *) mxCalloc(m+1,sizeof(int));
  iwsiz  = 7*m + 3;
  iwork  = (int *) mxCalloc(iwsiz,  sizeof(int));
  xlindx = (int *) mxCalloc(m+1,sizeof(int));
/* ------------------------------------------------------------
   Convert C-style symmetric matrix to adjacency structure
   (xadj,adjncy) in Fortran-style.
   ------------------------------------------------------------ */
  getadj(xadj,adjncy, Xjc,Xir,m);
/* ------------------------------------------------------------
   Convert PERM to integer, and make INVP
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++){
    j = permPr[i];
    perm[i] = j;
    invp[j-1] = i+1;
  }
/* ------------------------------------------------------------
   Initialize symbolic factorization
   Updates (PERM,INVP) to an equivalent ordering.
   ------------------------------------------------------------ */
  sfinit_(&m, Xjc+m,  xadj,adjncy, perm, invp, colcnt,
          &nnzl, &nsub, &nsuper, snode, xsuper, &iwsiz, iwork, &flag);
  if(flag == -1)
    mexErrMsgTxt("sfinit error.");
/* ------------------------------------------------------------
   Create output structure L
   ------------------------------------------------------------ */
  L_OUT = mxCreateStructMatrix(1, 1, NL_FIELDS, LFieldnames);
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); 
  outflag = mxGetPr(plhs[1]); 

/* ------------------------------------------------------------
   Create sparse output matrix L.L, m x m, with nnzl nonzeros.
   ------------------------------------------------------------ */
  L_FIELD = mxCreateSparse(m,m, nnzl,mxREAL);
  Ljc = mxGetJc(L_FIELD);
  Lir = mxGetIr(L_FIELD);
  Lpr = mxGetPr(L_FIELD);   
/* ------------------------------------------------------------
   Do symbolic factorization
   ------------------------------------------------------------ */
  symfct_(&m, Xjc+m,xadj,adjncy, perm,invp,colcnt,
          &nsuper,xsuper,snode, &nsub, xlindx, Lir, Ljc,
          &iwsiz,iwork, &flag);
  outflag[0] = (double)flag; 
  if(flag == -1)
    return;
    /*** mexErrMsgTxt("Insufficient working space."); ***/
  if(flag == -2)
    return;
    /*** mexErrMsgTxt("Input error symfct."); ***/
#ifdef DO_BFINIT
/* ------------------------------------------------------------
   Compute memory needs and cache-supernode-splitting for
   sparse block Cholesky
   ------------------------------------------------------------ */
  bfinit_(&m, &nsuper, xsuper,snode,xlindx, Lir,
          &cachsz, &tmpsiz, split);
#endif
/* ------------------------------------------------------------
   Expand row-indices from compact to standard subscript array,
   and fill nonzeros of L with 1's.
   ------------------------------------------------------------ */
  expandsub(m,nsuper,xsuper,xlindx,Ljc,Lir);
  for(i = 0; i < nnzl; i++)
    Lpr[i] = 1.0;
/* ------------------------------------------------------------
   Create output L.(L,perm,xsuper)
   ------------------------------------------------------------ */
  mxSetField(L_OUT, 0,"L", L_FIELD);                  /* L.L */
  L_FIELD = mxCreateDoubleMatrix(m, 1, mxREAL);       /* L.perm */
  permPr = mxGetPr(L_FIELD);
  mxSetField(L_OUT, 0,"perm", L_FIELD);
  L_FIELD = mxCreateDoubleMatrix(nsuper+1, 1, mxREAL);   /* L.xsuper */
  xsuperPr = mxGetPr(L_FIELD);
  mxSetField(L_OUT, 0,"xsuper", L_FIELD);
#ifdef DO_BFINIT
  L_FIELD = mxCreateDoubleMatrix(m, 1, mxREAL);          /* L.split */
  splitPr = mxGetPr(L_FIELD);
  mxSetField(L_OUT, 0,"split", L_FIELD);
  L_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL);          /* L.tmpsiz */
  *mxGetPr(L_FIELD) = tmpsiz;
  mxSetField(L_OUT, 0,"tmpsiz", L_FIELD);
#endif
/* ------------------------------------------------------------
   Convert (perm, xsuper) to floating point.
   ------------------------------------------------------------ */
  for(i = 0; i < m; i++)
    permPr[i] = perm[i];
  for(i = 0; i <= nsuper; i++)
    xsuperPr[i] = xsuper[i];
#ifdef DO_BFINIT
  for(i = 0; i < m; i++)
    splitPr[i] = split[i];
#endif
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(iwork);
  mxFree(invp);
  mxFree(perm);
  mxFree(xadj);
  mxFree(adjncy);
  mxFree(xlindx);
  mxFree(xsuper);
  mxFree(snode);
  mxFree(colcnt);
}
