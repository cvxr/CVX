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
#include "mex.h"
#include "blkchol.h"
/* ============================================================
   BACKWARD SOLVE:
   ============================================================ */
/* ************************************************************
   PROCEDURE bwsolve -- Solve y from L'*y = b, where
     L is lower-triangular.
   INPUT
     Ljc, Lir, Lpr - sparse lower triangular matrix
     xsuper - starting column in L for each (dense) supernode.
     nsuper - number of super nodes
   UPDATED
     y - full xsuper[nsuper]-vector, yOUTPUT = L' \ yINPUT.
   ************************************************************ */
void bwsolve(double *y, const mwIndex *Ljc, const mwIndex *Lir, const double *Lpr,
             const mwIndex *xsuper, const mwIndex nsuper)
{
  mwIndex jsup,j,inz,k;
  double yj, ljj;

  /* ------------------------------------------------------------
     For each supernode jsup:
     ------------------------------------------------------------ */
  j = xsuper[nsuper];      /* column after current snode (j=m)*/
  for(jsup = nsuper; jsup > 0; jsup--){
    for(k = 0; k < xsuper[jsup] - xsuper[jsup-1]; k++){
  /* ------------------------------------------------------------
     The equation L(:,j)'*yNEW = yOLD(j), yields
       y(j) -= L(j+1:m,j)'*y.
     ------------------------------------------------------------ */
      inz = Ljc[j-1];
      inz++;                        /* jump over diagonal entry */
      yj = realdot(Lpr+inz, y+j, k);
      for(inz += k; inz < Ljc[j]; inz++)
	yj += Lpr[inz] * y[Lir[inz]];
      y[--j] -= yj;
    }
  }
}

/* ************************************************************
   PROCEDURE partbwsolve -- Solve y from L(0:m-1,0:m-1)'*y = b, where
     L is lower-triangular.
   INPUT
     L=(ljc,lpr,xlindx,lindx) - sparse lower triangular matrix
     xsuper - starting column in L for each (dense) supernode.
     nsuper - number of super nodes (for NW-subblock)
     m  - order of L'- (sub) block, xsuper[nsuper-1] < m <= xsuper[nsuper]
       To solve with the complete L, choose m = xsuper[nsuper].
   UPDATED
     y - full m-vector, yOUT = L(0:m-1,0:m-1)'\yIN.
   ************************************************************ */
void partbwsolve(double *y, const mwIndex *ljc, const double *lpr,
             const mwIndex *xlindx, const mwIndex *lindx, const mwIndex *xsuper,
             const mwIndex nsuper, const mwIndex m)
{
  mwIndex jsup,j,inz,k,i, ixfirst,ixnz;
  double yj;
/* ------------------------------------------------------------
   For each supernode jsup:
   Let ixfirst point to the 1st row-subscript below current supernode.
   ------------------------------------------------------------ */
  j = m;           /* column/row after current entry */
  for(jsup = nsuper; jsup > 0; jsup--){
    ixfirst = xlindx[jsup-1] + (xsuper[jsup] - xsuper[jsup-1]);
/* ------------------------------------------------------------
   Case 1: L has sparse nonzeros below current supernode, but not
   beyond m:
   ------------------------------------------------------------ */
    if(ixfirst < xlindx[jsup])
      if(lindx[xlindx[jsup] - 1] < m)
        for(k = 0; j > xsuper[jsup-1]; k++){
  /* ------------------------------------------------------------
     The equation L(:,j)'*y=b(j), yields
       y(j) = b(j)-L(j+1:m,j)'*y.
     ------------------------------------------------------------ */
          inz = ljc[j-1] + 1;               /* jump over diagonal entry */
          yj = realdot(lpr+inz, y+j, k);
          for(inz += k, ixnz = ixfirst; inz < ljc[j]; inz++, ixnz++)
            yj += lpr[inz] * y[lindx[ixnz]];
          y[--j] -= yj;
        }
      else
/* ------------------------------------------------------------
   2nd CASE: L(:,j) has nonzeros beyond the requested order m,
   so use while i < m.
   ------------------------------------------------------------ */
        for(k = 0; j > xsuper[jsup-1]; k++){
          inz = ljc[j-1] + 1;               /* jump over diagonal entry */
          yj = realdot(lpr+inz, y+j, k);
          inz += k; ixnz = ixfirst;
          for(i=lindx[ixnz]; i < m; i = lindx[++ixnz])
            yj += lpr[inz++] * y[i];
          y[--j] -= yj;
        }
    else
/* ------------------------------------------------------------
   3rd CASE: L(:,j) has only nonzeros in snode jsup: purely dense.
   ------------------------------------------------------------ */
      for(k = 1, --j; j > xsuper[jsup-1]; k++){
        yj = realdot(lpr + ljc[j-1] + 1, y+j, k);
        y[--j] -= yj;
      }
  }
}
