/*%                                                          L = symbchol(X)
% SYMBCHOL Symbolic block sparse Cholesky factorization.
%   L = symbchol(X) returns a structure L that can be used
%   by the efficient block sparse Cholesky solver SPARCHOL.
%   The fields in L have the following meaning:
%
%   L.perm   - Multiple minimum degree ordering.
%
%   L.L      -  Sparse lower triangular matrix, has sparsity structure
%     of Cholesky factor of X(L.perm,L.perm).
%
%   L.xsuper - Supernode partition. Supernode jsup consists of
%     the nodes   L.xsuper(jsup) : L.xsuper(jsup)-1.
%
%   L.split  - Splitting of supernodes. Recommends to split supernode
%     in blocks of sizes   L.split(xsuper(jsup):L.xsuper(jsup)-1).
%
%   L.tmpsiz - Quantity used by SPARCHOL, to allocated enough working
%     storage.
%
%   L = symbchol(X,cachsz) optimizes L.split for a computer cache
%     of size CACHSZ * 1024 byte. Default cachsz = 16.
% SEE ALSO sparchol, fwblkslv, bwblkslv   

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

#define L_OUT    plhs[0]
#define NPAROUT 1

#define X_IN      prhs[0]
#define cachsz_IN prhs[1]
#define NPARINMIN 1
#define NPARIN 2


/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  mwIndex m, i,j, iwsiz, flag, nofsub, cachesiz, nsuper, nsub, nnzl;
  mwIndex *Xjc,*Xir, *snode, *xsuper, *invp, *colcnt, *xadj, *adjncy;
  mxArray *L_FIELD;
  const char *LFieldnames[] = {"L", "perm", "xsuper", "tmpsiz", "split"};
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  if(nrhs < NPARINMIN)
    mexErrMsgTxt("symbchol requires 1 or 2 input arguments");
  if(nlhs > NPAROUT)
    mexErrMsgTxt("symbchol produces 1 output argument");
/* ------------------------------------------------------------
   Check input X 
   ------------------------------------------------------------ */
  if(!mxIsSparse(X_IN))
    mexErrMsgTxt("Input matrix must be sparse");
  if( (m = mxGetM(X_IN)) != mxGetN(X_IN) )
    mexErrMsgTxt("X should be square.");
/* ------------------------------------------------------------
   Get input X
   ------------------------------------------------------------ */
  Xjc = mxGetJc(X_IN);
  Xir = mxGetIr(X_IN);
/* ------------------------------------------------------------
   Allocate working arrays:
   int xadj(m+1), adjncy(Xnnz), perm(m), invp(m), iwork(iwsiz)
   ------------------------------------------------------------ */
  xadj   = (mwIndex *) mxCalloc(m+1,sizeof(mwIndex));
  adjncy = (mwIndex *) mxCalloc(Xjc[m],sizeof(mwIndex));
  perm   = (mwIndex *) mxCalloc(m,sizeof(mwIndex));
  invp   = (mwIndex *) mxCalloc(m,sizeof(mwIndex));
  iwsiz  = 4 * m;
  iwork = (mwIndex *) mxCalloc(iwsiz,sizeof(mwIndex));
  
/* ------------------------------------------------------------
   Convert C-style symmetric matrix to adjacency structure
   (xadj,adjncy) in Fortran-style.
   ------------------------------------------------------------ */
  getadj(xadj,adjncy, Xjc,Xir,m);  
  /* ------------------------------------------------------------
   Compute multiple minimum degree ordering (J. Liu, in Fortran)
   ------------------------------------------------------------ */
  ordmmd_(&m,xadj,adjncy, invp,perm, &iwsiz,iwork, &nofsub, &flag);
  if(flag == -1)
    mexErrMsgTxt("Error in ordmmd.");
  
  mxFree(iwork);
  iwsiz  = 7*m + 3;
  iwork  = (mwIndex *) mxCalloc(iwsiz,  sizeof(mwIndex));
  colcnt = (mwIndex *) mxCalloc(m,  sizeof(mwIndex));
  snode  = (mwIndex *) mxCalloc(m,  sizeof(mwIndex));
  xsuper = (mwIndex *) mxCalloc(m+1,sizeof(mwIndex));
  xlindx = (mwIndex *) mxCalloc(m+1,sizeof(mwIndex));
  
  sfinit_(&m, Xjc+m,  xadj,adjncy, perm, invp, colcnt,
          &nnzl, &nsub, &nsuper, snode, xsuper, &iwsiz, iwork, &flag);
  if(flag == -1)
    mexErrMsgTxt("sfinit error.");
  
  /* ------------------------------------------------------------
   Do symbolic factorization
   ------------------------------------------------------------ */
  symfct_(&m, Xjc+m,xadj,adjncy, perm,invp,colcnt,
          &nsuper,xsuper,snode, &nsub, xlindx, Lir, Ljc,
          &iwsiz,iwork, &flag);
  if(flag == -1)
    mexErrMsgTxt("Insufficient working space.");
  if(flag == -2)
    mexErrMsgTxt("Input error symfct.");
  
  #ifdef DO_BFINIT
/* ------------------------------------------------------------
   Compute memory needs and cache-supernode-splitting for
   sparse block Cholesky
   ------------------------------------------------------------ */
  bfinit_(&m, &nsuper, xsuper,snode,xlindx, Lir,
          &cachsz, &tmpsiz, split);
#endif


  
}

/* ============================================================
   SUBROUTINES:
   ============================================================ */
/* ------------------------------------------------------------
   GETADJ - Copies off-diagonal entries from C-style sparse
     matrix (cjc,cir) to Fortran style sparse matrix (forjc,forir).
     On input, n is number of columns.
   ------------------------------------------------------------ */
void getadj(mwIndex *forjc,mwIndex *forir,const mwIndex *cjc,const mwIndex *cir,const mwIndex n)
{
	mwIndex i,j,inz,ix;

	inz = 0;
    for(j = 0; j < n; j++){
		forjc[j] = inz + 1;
		for(ix = cjc[j]; ix < cjc[j+1]; ix++)
			if((i = cir[ix]) != j)
				forir[inz++] = ++i;
	}
	forjc[n] = ++inz;
}
