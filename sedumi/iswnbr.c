/*
%                                 [delta,h,alpha] = iswnbr(vSQR,thetaSQR)
% ISWNBR  Checks feasibility w.r.t. wide region/neighborhood of Sturm-Zhang.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function [delta,h,alpha] = iswnbr(vSQR,thetaSQR)

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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"

#define	DELTA_OUT myplhs[0]
#define	H_OUT myplhs[1]
#define	ALPHA_OUT myplhs[2]
#define NPAROUT 3

#define	VSQR_IN prhs[0]
#define	THETASQR_IN prhs[1]
#define NPARIN 2

/* ------------------------------------------------------------
   Macros:
   ------------------------------------------------------------ */
#if !defined(SQR)
#define SQR(x) ((x)*(x))
#endif

#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

/* --------------------------------------
   FLOAT COMPARE: FOR SORTING A FLOAT ARRAY
   -------------------------------------- */
typedef int (*COMPFUN)(const void *pa,const void *pb);
#define fsort(vec,n)  qsort((void *)(vec), (n), sizeof(double), (COMPFUN) fcmp);

signed char fcmp(const double *a, const double *b)
{
   return( (*a > *b) - (*a < *b)  );
}

/********   GETDELTA  .... THE WIDE NBRHD MEMBERSHIP TEST   ********
     PURPOSE - THIS ROUTINE COMPUTES THE PROXIMITY MEASURE
	WITH RESPECT TO THE WIDE REGION C(THETA). IT KEEPS A
	GROWING SUBSET AND A SHRINKING SUPERSET OF T IN A
	LINKED LIST, SO THAT THE ITERATIVE CALCULATIONS
	ARE MINIMIZED.
     INPUT:
       vSQR - full n x 1 vector of squared v-space solution.
       thetaSQR - Squared parameter of region C(theta).
       n - order, i.e. length(vSQR).
     OUTPUT:
       h, alpha - such that vTAR = (1-alpha)*max(h,v) is the
          projection of v onto theta-central region.
     RETURNS:
       delta - sqrt(r) * norm(vTAR - v) / norm(v),   with r := n / theta^2 
         (i.e., sine measure)
         Returns 1e100 if w(j) <= 0 for some j.
     WORKING ARRAY
       wQ - length n vector of doubles.
 ****************************************************************/

double getdelta(double *ph, double *palpha, const double *w,
		const double thetaSQR, const mwIndex n, double *wQ)
{
 double gap,r,h,hSQR,oldhSQR,hubSQR, sumdifv,sumdifw, sumwNT,wj,
   deltaSQR,alpha;
 mwIndex cardT,cardQ,i,j,STOP;

 /* ------------------------------------------------------------
    gap = sum(w),  r = n / theta^2
    ------------------------------------------------------------ */
 for(i = 0, gap = 0.0; i < n; i++)
   gap += w[i];
 r = n / thetaSQR;
 /* ------------------------------------------------------------
    In the following, T is the index set of w[i]'s that are too small,
    and Q is the set for which we don't know in the first pass.
    However, after sorting wQ=w[Q], we will sort that out too.
    ------------------------------------------------------------ */
/* ------------------------------------------------------------------
   IF THETA == 1, WE HAVE:
     h = sqrt(max(w))
     sumdifv = sum (h-v_i) = n*h - sum(sqrt(w)),   (v := sqrt(w))
     sumdifw = sum (h^2-w_i) = n*h^2 - gap.
   Use 2 passes to compute sumdifv,sumdifw in stable way.
   ------------------------------------------------------------------ */
 if(1.0 - thetaSQR <= 1E-8){
   for(i = 0, hSQR = 0.0; i < n; i++)
     hSQR = MAX(hSQR,w[i]);
   h = sqrt(hSQR);
   sumdifv = 0.0; sumdifw = 0.0;
   for(i = 0; i < n; i++){
     sumdifw += hSQR - w[i];
     sumdifv += h - sqrt(w[i]);
   }
 }
 else{
/* ------------------------------------------------------------
   0 < THETA < 1:
     LB: hSQR = sumwNT/(r-|T|)
     UB: sumwNT/(r-|T| - |Q|)         (Q is first stage inconclusive set)
    sumdifv = sum_{j in T} (h-v_j)
    sumdifw = sum_{j in T} (h-w_j)
  Notice that sum(dif) is much stabler than dif(sum).
  i : next entry for wQ
  ------------------------------------------------------------ */
   sumwNT = gap;
   cardT = 0; sumdifv = 0; sumdifw = 0;
   cardQ = n;
   i = 0;
   hSQR = sumwNT / (r - cardT);
   hubSQR = sumwNT / (r-(n-1));
   for(j = 0; j < n; j++){
     wj = w[j];
     if(wj >= hubSQR){               /* wj >= hubSQR ==> not in T */
       --cardQ;
       hubSQR = sumwNT / (r-cardT-cardQ);
     }
     else if(wj < hSQR){             /* wj < hSQR ==> in T */
       if(wj <= 0.0)
         return 1e100;                 /* error: w should be positive */
       ++cardT;
       mxAssert(cardQ>0,"");
       --cardQ;
       hubSQR *= 1 - wj/sumwNT;
       sumwNT -= wj;
       oldhSQR = hSQR;
       hSQR = sumwNT / (r - cardT);
       sumdifw += (oldhSQR-wj) + cardT * (hSQR-oldhSQR);
       sumdifv += (sqrt(oldhSQR)-sqrt(wj))+ cardT*(sqrt(hSQR)-sqrt(oldhSQR));
     }
     else                            /* Inconclusive: j in Q */
       wQ[i++] = wj;
   }
   mxAssert(i == cardQ,"");
   /* ------------------------------------------------------------
      The same treatment for the Q set, but we
      sort the (presumably short) wQ first.
     ------------------------------------------------------------ */
   if(cardQ){
     fsort(wQ, cardQ);
     for(STOP = 0, j = 0; !STOP; ){
       wj = wQ[j];
       if(wj >= hSQR)
	 STOP = 1;
       else{
	 ++cardT;
	 sumwNT -= wj;
	 oldhSQR = hSQR;
	 hSQR = sumwNT / (r - cardT);
	 sumdifw += (oldhSQR-wj) + cardT * (hSQR-oldhSQR);
	 sumdifv += (sqrt(oldhSQR)-sqrt(wj)) +
	   cardT * (sqrt(hSQR)-sqrt(oldhSQR));
	 STOP = (cardQ == ++j);
       }
     }
   } /* cardQ > 0 */
   /* ------------------------------------------------------------
      Let h := sqrt(hSQR)
      ------------------------------------------------------------ */
   h = sqrt(hSQR);
 }  /* theta != 1 */
 /* ------------------------------------------------------------
    FOR ALL THETA :
      alpha = sumdifv/(r*h)
      deltaSQR = r * ( 2*alpha-alpha^2 - (1-alpha)^2 * sumdifw/gap )
    (THE ABOVE DIFFERENCE SHOULD NOT BE NUMERICALLY DANGEROUS,
     SINCE alpha IS *SIGNIF* BIGGER THAN sumdifw/gap )
    ------------------------------------------------------------ */
 alpha = sumdifv/ (r*h);
 deltaSQR = alpha*(2-alpha) - SQR(1-alpha) * sumdifw/gap;
 *palpha = alpha;
 *ph = h;
 return sqrt(r * deltaSQR);
}

/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
     [delta,h,alpha] = iswnbr(vSQR,thetaSQR)
     Computes proximity "delta" w.r.t cregion C(theta).
     The projection of v onto C(theta) is (1-alpha)*max(h,v).
   ************************************************************ */
void mexFunction(int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
 mxArray *myplhs[NPAROUT];
 double *w,*pdelta,*ph,*palpha, *fwork;
 double thetaSQR;
 mwIndex i,n;

/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "iswnbr requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "iswnbr produces less output arguments");
 /* ------------------------------------------------------------
    Get input vector w:=vSQR and cregion parameter thetaSQR
    ------------------------------------------------------------ */
 w = mxGetPr(VSQR_IN);
 n = mxGetM(VSQR_IN) * mxGetN(VSQR_IN);
 thetaSQR = mxGetScalar(THETASQR_IN);
 /* ------------------------------------------------------------
    Allocate output DELTA, H, ALPHA.
    ------------------------------------------------------------ */
 DELTA_OUT = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
 pdelta = mxGetPr(DELTA_OUT);
 H_OUT = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
 ph = mxGetPr(H_OUT);
 ALPHA_OUT = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
 palpha = mxGetPr(ALPHA_OUT);
 /* ------------------------------------------------------------
    Allocate working array fwork(n)
    ------------------------------------------------------------ */
 fwork = (double *) mxCalloc(n,sizeof(double));
 /* ------------------------------------------------------------
    The actual job is done here:.
    ------------------------------------------------------------ */
 *pdelta = getdelta(ph,palpha, w,thetaSQR,n, fwork);
 /* ------------------------------------------------------------
    RELEASE WORKING ARRAY.
    ------------------------------------------------------------ */
 mxFree(fwork);
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  i = MAX(nlhs, 1);
  memcpy(plhs,myplhs, i * sizeof(mxArray *));
  for(; i < NPAROUT; i++)
    mxDestroyArray(myplhs[i]);
}
