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


/* ============================================================
   DOT-PRODUCT AND VECTOR ARRAY-OPS.
   ============================================================ */

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- realHadamard
   Computes  r = x .* y  using loop-unrolling.
   ************************************************************ */
void realHadamard(double * r, const double *x, const double *y, const mwIndex n)
{
/* TODO: This blas solution segfaults. OpenMP?
//     int one=1;
//     double zero=0.0;
//     #ifdef PC
//     dcopy(&n,y,&one,r,&one);
//     dtbmv('u','n','n',&n,&zero,x,&one,r,&one);
//     #endif
//     
//     #ifdef UNIX
//     dcopy_(&n,y,&one,r,&one);
//     dtbmv_('u','n','n',&n,&zero,x,&one,r,&one);
//     #endif
//     return;*/
    mwIndex i;

 for(i=0; i+7< n; ){              /* LEVEL 8 */
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
 }
 if(i+3 < n){                         /* LEVEL 4 */
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
 }
/* ------------------------------------------------------------
   Now, i in {n-3, n-2, n-1, n}. Do the last n-i elements.
   ------------------------------------------------------------ */
 if(i+1 < n){                        /* LEVEL 2 */
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
 }
 if(i< n)                            /* LEVEL 1 */
   r[i] = x[i] * y[i];
}

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- realHadadiv
   Computes  r = x ./ y  using loop-unrolling.
   ************************************************************ */
void realHadadiv(double * r, const double *x, const double *y, const mwIndex n)
{
 mwIndex i;

 for(i=0; i+3< n; ){              /* LEVEL 4 */
   r[i] = x[i] / y[i]; i++;
   r[i] = x[i] / y[i]; i++;
   r[i] = x[i] / y[i]; i++;
   r[i] = x[i] / y[i]; i++;
 }
/* ------------------------------------------------------------
   Now, i in {n-3, n-2, n-1, n}. Do the last n-i elements.
   ------------------------------------------------------------ */
 if(i+1 < n){                        /* LEVEL 2 */
   r[i] = x[i] / y[i]; i++;
   r[i] = x[i] / y[i]; i++;
 }
 if(i< n)                            /* LEVEL 1 */
   r[i] = x[i] / y[i];
}

