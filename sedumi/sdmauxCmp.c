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

#include <string.h>     /* memcpy etc */
#include <math.h>      /* floor and log */
#include "blksdp.h"

/* ************************************************************
   COMPARISON/ SEARCHING
   ************************************************************ */
/* ------------------------------------------------------------
   For ascending sort (e.g. 0,1,2,.. or a,b,c,..) : 
   compare yields < 0 iff a < b, 0 iff a==b, > 0 iff a > b.
   ------------------------------------------------------------ */
/* Integer compare: (for ibsearch) */
char icmp(const mwIndex *a, const mwIndex *b)
{
   return( (*a > *b) - (*a < *b)  );
}

/* Integer compare (with integer key): (for kiqsort) */
char kicmp(const keyint *a, const keyint *b)
{
   return( (a->i > b->i) - (a->i < b->i)  );
}

/* Float compare (with integer key), for DESCENDING sort (e.g. 10,9,8,..).*/
char kdcmpdec(const keydouble *a, const keydouble *b)
{
   return( (a->r < b->r) - (a->r > b->r)  );
}

/* ************************************************************
   PROCEDURE intbsearch - search through (in strictly ascending order)
      sorted integer array.  DIFFERENT FROM ansi BSEARCH, since it returns
      always the offset k=*pi, such that x(0:k-1) < key <= x(k:n-1).
      Thus, x(k) CAN be larger than key (viz. if key not in array).
      if x(0:n-1) < key then k=n.
   INPUT
     x - length n integer array
     key - integer to search for
     n - length of x
   UPDATED
     pi - On input, 0<= *pi; we'll search only in x[*pi,n-1].
       On output, x(0:*pi-1) < key, x(*pi:n) >= key.
       NB1: *piIN <= *piOUT <= n, so possibly *piOUT=n.
       NB2: if *piIN==n then *piOUT=n.
   RETURNS
     1 if found, 0 otherwise. If found, then x[*pi]=key, otherwise
     x[*pi] > key.
   ************************************************************ */
bool intbsearch(mwIndex *pi, const mwIndex *x, const mwIndex n, const mwIndex key)
{
 mwIndex i,j,r;

 i = *pi;
 mxAssert(i >= 0,"");
 if(i < n){
   if(x[i] < key){
     r = n;
/* ------------------------------------------------------------
   During the loop, x[i] < key and r<n => key < x[r], i.e. key has
   to be strictly
   within (i,r). Therefore, we also take j strictly in (i,r).
   ------------------------------------------------------------ */
     for(j = (i+n) / 2; j > i; j = (i+r) / 2){
       if(x[j] > key)
         r = j;                    /* new right limit */
       else if(x[j] < key)
         i = j;                    /* new left limit */
       else{
         *pi = j;                  /* Found : x[j] == key */
         return 1;
       }
     }
     *pi = r;               /* x[r] > key */
     return 0;              /* not found: i==j==(r-1), x[r-1] < key < x[r] */
   }
/* ------------------------------------------------------------
   If, for the initial i, x[i] >= key, then keep this i, and
   return found = (x[i] == key).
   ------------------------------------------------------------ */
   else /* x[(i = *pi)] >= key */
     return (x[i] == key);
 }
 return 0;     /* if i>=n, i.e. the list is empty */
}

/* ************************************************************
   PROCEDURE intmbsearch - multi-bsearch: locate positions of
     ascending array y in the ascending array x. If y is a singleton,
     then this corresponds to regular bsearch.
   INPUT
     x, y - integer arrays of length xnnz and ynnz, resp.
     xnnz, ynnz - lengths of x and y, resp.
     bsize - size of mwIndex working array b; bsize = floor(log(ynnz+1)/log(2)).
   OUTPUT
     z - length ynnz+2 array. Sets z[0]=0, z[ynnz+1]=xnnz, and
       x[z[k+1]-1] < y[k] <= x[z[k+1]] for all k=0:ynnz-1.
     found - length ynnz array. found[k] = (x[z[k+1]] == y[k]).
   WORKING ARRAY
     b - length bsize iwork array; bsize = floor(log(ynnz+1)/log(2)).
   RETURNS 0 SUCCESS, -1 bsize too small.
   ************************************************************ */
char intmbsearch(mwIndex *z, bool *found, const mwIndex *x, const mwIndex xnnz,
		const mwIndex *y, const mwIndex ynnz, mwIndex *b, const mwIndex bsize)
{
  mwIndex a,bk,k,m;
  mwIndex *zp1;
  bool isknegative;
/* ------------------------------------------------------------
   Init:
   ------------------------------------------------------------ */
  mxAssert(xnnz >= 0 && ynnz >= 0,"");
  if(bsize < (mwIndex) floor(log(ynnz+1.0) / log(2.0)))
    return -1;                   /* ERROR: insufficient work space */
  z[0] = 0;
  zp1 = z + 1;
  zp1[ynnz] = xnnz;
  if(ynnz == 0)
    return 0;
/* ------------------------------------------------------------
   Divide and Conquer:
   Divide y in binary way; k is current depth in tree;
   Still need to locate [y[a], y[bk]); the uncertainty interval
   in x is [z[a], zp1[bk]).
   Note that we located x[zp1[bk]-1] < y[bk] <= x[zp1[bk]].
   All y[0:a-1] have already been located: x[z[a]-1] < y[a-1].
   ------------------------------------------------------------ */
  k = 0;
  a = 0;
  b[0] = ynnz;
  isknegative=0;
  while(!isknegative){
    bk = b[k];
    mxAssert(a < bk,"");
/* ------------------------------------------------------------
   y[a:bk-1] have not yet been located. Locate y[m] with
   a <= m < bk
   ------------------------------------------------------------ */
    m = (mwIndex) floor((double) (a+bk-1)/2); 
    zp1[m] = z[a];                        /* must be in [z[a],zp1[bk]) */
    found[m] = intbsearch(zp1+m, x, zp1[bk], y[m]);
/* ------------------------------------------------------------
   If a < m then m < bk-1 and we DESCEND in the tree to locate y[a,m-1].
   ------------------------------------------------------------ */
    if(m > a){
      b[++k] = m;                        /* Locate y[a:m-1] */
      mxAssert(k < bsize,"");            /* Also b[k] < bk - 1 */
    }
/* ------------------------------------------------------------
   If a = m  = bk-2, then STAY on this level with a = m+1, i.e. ++a.
   We've located y[0:m].
   ------------------------------------------------------------ */
    else if(++a >= bk){
/* ------------------------------------------------------------
   If a:bk-1 = {a}, i.e. m=a, then we finished this branch of the tree.
   We've located y[0:bk]. ASCEND in tree with a = bk+1.
   ------------------------------------------------------------ */
      if(k>0)
        --k;
      else
          isknegative=1;
      mxAssert(k>=0,"");
      ++a;         /* Note that a = bk+1 < b[k], see DESCEND */
      mxAssert(a == bk+1,"");
    }
  } /* if k >=0 then proceed location on y[a,b[k]-1]. */
/* ------------------------------------------------------------
   SUCCESS: return 0.
   ------------------------------------------------------------ */
  return 0;
}
