/* ************************************************************
   MODULE sdmaux*.c  -- Several low-level subroutines for the
   mex-files in the Self-Dual-Minimization package.

    This file is part of SeDuMi 1.04BETA
    Copyright (C) 2000 Jos F. Sturm
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
   ************************************************************ */

#include "blksdp.h"

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- r=realdot(x,y,n)
   Computes r=sum(x_i * y_i) using loop-unrolling.
   ************************************************************ */
double realdot(const double *x, const double *y, const int n)
{
 int i;
 double r;

 r = 0.0;
 for(i = 0; i < n-7; i++){          /* LEVEL 8 */
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i];
 }
/* ------------------------------------------------------------
   Now, i in {n-7, n-6, ..., n}. Do the last n-i elements.
   ------------------------------------------------------------ */
 if(i < n-3){                            /* LEVEL 4 */
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i]; i++;
 }
 if(i < n-1){                           /* LEVEL 2 */
   r+= x[i] * y[i]; i++;
   r+= x[i] * y[i]; i++;
 }
 if(i < n)                              /* LEVEL 1 */
   r+= x[i] * y[i];
 return r;
}

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- r=realssqr(x,n)
   Computes r=sum(x_i^2) using loop-unrolling.
   ************************************************************ */
double realssqr(const double *x, const int n)
{
  int i;
  double r;

  r=0.0;
  for(r=0.0, i=0; i< n-7; i++){          /* LEVEL 8 */
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]);
  }
/* ------------------------------------------------------------
   Now, i in {n-7, n-6, ..., n}. Do the last n-i elements.
   ------------------------------------------------------------ */
  if(i < n-3){                              /* LEVEL 4 */
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]); i++;
  }
  if(i < n-1){                           /* LEVEL 2 */
    r+= SQR(x[i]); i++;
    r+= SQR(x[i]); i++;
  }
  if(i < n)                              /* LEVEL 1 */
    r+= SQR(x[i]);
  return r;
}
