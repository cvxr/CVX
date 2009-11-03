/* ************************************************************
   MODULE sdmaux*.c  -- Several low-level subroutines for the
   mex-files in the Self-Dual-Minimization package.

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

#include "blksdp.h"
/* ************************************************************
   TIME-CRITICAL PROCEDURE -- scalarmul
   Computes  r = alpha * x  using BLAS.
   ************************************************************ */
void scalarmul(double *r, const double alpha,const double *x,const mwIndex n)
{
  mwIndex one=1;
  #ifdef PC
  dcopy(&n,x,&one,r,&one);
  dscal(&n,&alpha,r,&one);
  #endif
  #ifdef UNIX
  dcopy_(&n,x,&one,r,&one);
  dscal_(&n,&alpha,r,&one);
  #endif  
  return;
  
}

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- addscalarmul
   Computes  r += alpha * x  using BLAS.
   ************************************************************ */
void addscalarmul(double *r, const double alpha,const double *x,const mwIndex n)
{
  /*USE BLAS*/
    mwIndex one=1;
    #ifdef PC
    daxpy(&n,&alpha,x,&one,r,&one);
    #endif
    #ifdef UNIX
    daxpy_(&n,&alpha,x,&one,r,&one);
    #endif    
    return;
  
}

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- subscalarmul(x,alpha,y,n)
   Computes x -= alpha * y using BLAS.
   ************************************************************ */
void subscalarmul(double *x, const double alpha, const double *y, const mwIndex n)
{
  /*USE BLAS*/
    mwIndex one=1;
    const double minusalpha=-alpha;
    #ifdef PC
    daxpy(&n,&minusalpha,y,&one,x,&one);
    #endif
    #ifdef UNIX
    daxpy_(&n,&minusalpha,y,&one,x,&one);
    #endif    
    return;
}
