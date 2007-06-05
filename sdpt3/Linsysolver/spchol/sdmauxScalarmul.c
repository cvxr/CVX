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
   TIME-CRITICAL PROCEDURE -- scalarmul
   Computes  r = alpha * x  using loop-unrolling.
   ************************************************************ */
void scalarmul(double *r, const double alpha,const double *x,const int n)
{
  int k;

  for(k = 0; k < n-7; k++){                 /* LEVEL 8 */
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k];
  }
/* ------------------------------------------------------------
   Now, i in {n-7, n-6, ..., n}. Do the last n-i elements.
   ------------------------------------------------------------ */
  if(k < n-3){                              /* LEVEL 4 */
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k]; k++;
  }
  if(k < n-1){                              /* LEVEL 2 */
    r[k] = alpha * x[k]; k++;
    r[k] = alpha * x[k]; k++;
  }
  if(k < n)                                 /* LEVEL 1 */
    r[k] = alpha * x[k];
}

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- addscalarmul
   Computes  r += alpha * x  using loop-unrolling.
   ************************************************************ */
void addscalarmul(double *r, const double alpha,const double *x,const int n)
{
  int k;

  for(k = 0; k < n-7; k++){                 /* LEVEL 4 */
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k];
  }
/* ------------------------------------------------------------
   Now, i in {n-7, n-6, ..., n}. Do the last n-i elements.
   ------------------------------------------------------------ */
  if(k < n-3){                              /* LEVEL 4 */
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k]; k++;
  }
  if(k < n-1){                              /* LEVEL 2 */
    r[k] += alpha * x[k]; k++;
    r[k] += alpha * x[k]; k++;
  }
  if(k < n)                                 /* LEVEL 1 */
    r[k] += alpha * x[k];
}

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- subscalarmul(x,alpha,y,n)
   Computes x -= alpha * y using LEVEL 8 loop-unrolling.
   ************************************************************ */
void subscalarmul(double *x, const double alpha, const double *y, const int n)
{
  int i;
  
  for(i=0; i< n-7; i++){          /* LEVEL 8 */
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i];
  }
/* ------------------------------------------------------------
   Now, i in {n-7, n-6, ..., n}. Do the last n-i elements.
   ------------------------------------------------------------ */
  if(i < n-3){                           /* LEVEL 4 */
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
  }
  if(i < n-1){                           /* LEVEL 2 */
    x[i] -= alpha * y[i]; i++;
    x[i] -= alpha * y[i]; i++;
  }
  if(i < n)                              /* LEVEL 1 */
    x[i] -= alpha * y[i];
}
