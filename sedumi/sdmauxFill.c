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
/* ============================================================
   ARRAY FILLING
   ============================================================ */
/* ************************************************************
   TIME CRITICAL PROCEDURE fromto -- x = i:n-1
   INPUT
     i, n - first and just beyond last integer to write in x.
   OUTPUT
     x - (n-i) array, will contain x = i:(n-1)
   ************************************************************ */
void fromto(mwIndex *x, mwIndex i, const mwIndex n)
{
  x -= i;                /* makes x(i:n-1) valid */
  for(; i < n; i++)
    x[i] = i;
}

/* ************************************************************
   PROCEDURE fzeros -- z = zeros(n,1)
   INPUT  n = length(z)
   OUTPUT z = zeros(n,1)
   ************************************************************ */
void fzeros(double *z,const mwIndex n)
{
  mwIndex k;
  
  for(k=0; k+3 < n; ){           /* LEVEL 4 */
    z[k] = 0.0; k++;
    z[k] = 0.0; k++;
    z[k] = 0.0; k++;
    z[k] = 0.0; k++;
  }
  if(k+1 < n){                      /* LEVEL 2 */
    z[k] = 0.0; k++;
    z[k] = 0.0; k++;
  }
  if(k < n)                           /* LEVEL 1 */
    z[k] = 0.0;
}
