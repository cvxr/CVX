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

#include <string.h>
#include "blksdp.h"
/* ************************************************************
   PROCEDURE: triudotprod - Computes trace(X*Y), where X and Y are
     supposed to be symmetric; only the triu's are used.
   INPUT
     x - full n x n, only triu(x) used
     y - full n x n, only triu(x) used
     n - order
   RETURNS the inner product between X and Y, trace(X*Y).
   ************************************************************ */
double triudotprod(const double *x, const double *y, const mwIndex n)
{
  mwIndex i, j;
  double z;
  if(n <= 0)
    return 0.0;
  z = x[0] * y[0];
  for(j = 1; j < n; j++){
    x += n; y += n;                  /* point to x(:,j) and y(:,j) */
    z += 2 * realdot(x,y,j) + x[j] * y[j];
  }
  return z;
}

/* striudotprod: z = tr(X*Y), where diag(X)=0. only triu(X,1) and triu(Y,1)
   are used. For the imaginary part into the real inner product. */
double striudotprod(const double *x, const double *y, const mwIndex n)
{
  mwIndex i, j;
  double z;
  if(n <= 1)
    return 0.0;
  z = 0.0;
  for(j = 1; j < n; j++){
    x += n; y += n;                  /* point to x(:,j) and y(:,j) */
    z += realdot(x,y,j);
  }
  return 2 * z;
}
