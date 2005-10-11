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
void realHadamard(double * r, const double *x, const double *y, const int n)
{
 int i;

 for(i=0; i< n-7; ){              /* LEVEL 8 */
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
 }
 if(i < n-3){                         /* LEVEL 4 */
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
   r[i] = x[i] * y[i]; i++;
 }
/* ------------------------------------------------------------
   Now, i in {n-3, n-2, n-1, n}. Do the last n-i elements.
   ------------------------------------------------------------ */
 if(i < n-1){                        /* LEVEL 2 */
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
void realHadadiv(double * r, const double *x, const double *y, const int n)
{
 int i;

 for(i=0; i< n-3; ){              /* LEVEL 4 */
   r[i] = x[i] / y[i]; i++;
   r[i] = x[i] / y[i]; i++;
   r[i] = x[i] / y[i]; i++;
   r[i] = x[i] / y[i]; i++;
 }
/* ------------------------------------------------------------
   Now, i in {n-3, n-2, n-1, n}. Do the last n-i elements.
   ------------------------------------------------------------ */
 if(i < n-1){                        /* LEVEL 2 */
   r[i] = x[i] / y[i]; i++;
   r[i] = x[i] / y[i]; i++;
 }
 if(i< n)                            /* LEVEL 1 */
   r[i] = x[i] / y[i];
}

#ifdef SEDUMI_OLD
/* ============================================================
   LORENTZ OPERATIONS
   ============================================================ */
/* ************************************************************
   PROCEDURE qscale : LORENTZ SCALE D(x)y = z + mu * x (full version)
     mu = (y1+alpha)/sqrt(2), z = rdetx * [alpha; y(2:n)],
     where alpha = (x(2:n)'*y(2:n)) / (x(1)+ sqrt(2) * rdetx)
   INPUT
     x,y - full n x 1
     rdetx - sqrt(det(x))
     n - order of x,y,z.
   OUTPUT
     z - full n x 1. Let z := rdetx * [alpha; y(2:n)].
   RETURNS
     mu = (y1+alpha)/sqrt(2).
   ************************************************************ */
double qscale(double *z,const double *x,const double *y,
              const double rdetx,const int n)
{
 double alpha, mu;
/* ------------------------------------------------------------
   alpha = (x(2:n)'*y(2:n)) / (x(1)+ sqrt(2) * rdetx)
   ------------------------------------------------------------ */
 alpha = realdot(x+1,y+1,n-1) / (x[0] + M_SQRT2 * rdetx);
/* ------------------------------------------------------------
   z = rdetx * [alpha; y(2:n)].
   ------------------------------------------------------------ */
 z[0] = rdetx * alpha;
 scalarmul(z+1,rdetx,y+1,n-1);
/* ------------------------------------------------------------
   RETURN mu = (y1+alpha)/sqrt(2).
   ------------------------------------------------------------ */
 return (y[0] + alpha) / M_SQRT2;
}

/* ************************************************************
   PROCEDURE qlmul : LORENTZ SCALE z = D(x)y (full version)
     z=D(x)y = [x'*y / sqrt(2);  mu * x(2:n) + rdetx * y(2:n)],
     where mu = (z(1)+rdetx*y1) / (x(1)+ sqrt(2) * rdetx)
   INPUT
     x,y - full n x 1
     rdetx - sqrt(det(x))
     n - order of x,y,z.
   OUTPUT
     z - full n x 1. Let z := D(x)y.
   ************************************************************ */
void qlmul(double *z,const double *x,const double *y,
	   const double rdetx,const int n)
{
  double mu;
  mu = qscale(z, x,y,rdetx,n);
  addscalarmul(z,mu,x,n);
}

/* ************************************************************
   PROCEDURE qldiv : LORENTZ SCALE z = D(x)\y (full version)
    D(x)\y = (1/det x) * [x'Jy/sqrt(2); rdetx * y2-alpha*x2],
    where alpha = (x'Jy/sqrt(2) + rdetx*y1) / (x(1)+ sqrt(2) * rdetx)
   INPUT
     x,y - full n x 1
     rdetx - sqrt(det(x))
     n - order of x,y,z.
   OUTPUT
     z - full n x 1. Let z := D(x)^{-1}y.
   ************************************************************ */
void qldiv(double *z,const double *x,const double *y,
	   const double rdetx,const int n)
{
 double alpha,x1,y1,z1;
/* ------------------------------------------------------------
   z1 = x'*J*y / (sqrt(2) * det x),
   alpha = (z1+y1/rdetx) / (x(1)+ sqrt(2) * rdetx)
   ------------------------------------------------------------ */
 x1 = x[0]; y1 = y[0];
 z1 = (x1*y1 - realdot(x+1,y+1,n-1)) / (M_SQRT2 * SQR(rdetx));
 alpha = (z1 + y1 / rdetx) / (x1 + M_SQRT2 * rdetx);
/* ------------------------------------------------------------
   z(1) = z1, z(2:n) = y(2:n)/rdetx - alpha * x(2:n).
   ------------------------------------------------------------ */
 z[0] = z1;
 scalarmul(z+1,-alpha,x+1,n-1);
 addscalarmul(z+1,1/rdetx,y+1,n-1);
}
#endif
