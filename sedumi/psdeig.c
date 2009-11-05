/*
%                                                [lab,q] = psdeig(x,K)
% PSDEIG  Computes spectral coefficients of x w.r.t. K
%   Arguments "q" is optional - without it's considerably faster.
%   FLOPS indication: 1.3 nk^3 versus 9.0 nk^3 for nk=500,
%                     1.5 nk^3        9.8 nk^3 for nk=50.
%
% SEE ALSO sedumi
% **********  INTERNAL FUNCTION OF SEDUMI **********
function [lab,q] = psdeig(x,K)
 
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
#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define USE_SVD 1 /* Uses [q,lab,q]=svd(X) for complex hermitian (psd?) x */

#define LAB_OUT plhs[0]
#define Q_OUT plhs[1]
#define NPAROUT 2

#define X_IN prhs[0]
#define K_IN prhs[1]
#define NPARIN 2

/*
 * SYMPROJ: projection onto the symmetric matrices
 * INPUT:   x, n: a full n x n matrix x
 * OUTPUT:  y: contains (x+x')/2
 */

void symproj( double *y, const double *x, const mwSize n )
{
    mwIndex colp, i, j;
    double yij;
    for ( i = 0 ; i != n ; x += n, y += n, ++i ) {
        y[i] = x[i];
        for ( colp = n + i, j = i + 1 ; j != n ; ++j, colp += n ) {
            yij = ( x[j] + x[colp] ) / 2;
            y[j] = yij;
            y[colp] = yij;
        }
    }
}

/*
 * SKEWPROJ: projection onto skew-symmetric matrices
 * INPUT:    x, n: full n x n matrix x
 * OUTPUT:   y: contains (x-x')/2
 */

void skewproj( double *y, const double *x, const mwSize n )
{
    mwIndex colp,i,j;
    double yij;
    for ( i = 0 ; i != n ; x += n, y += n, ++i ) {
        y[i] = 0;
        for ( colp = n + i, j = i + 1 ; j != n ; ++j, colp += n ) {
            yij = ( x[j] - x[colp] ) / 2;
            y[j] = yij;
            y[colp] = -yij;
        }
    }
}

/*
 * MEXFUNCTION: main Matlab entry point
 */

void mexFunction(
    const int nlhs, mxArray *plhs[],
    const int nrhs, const mxArray *prhs[] )
{
    mxArray *output_array[3], *Xk;
    coneK cK;
    mwSize nk, nkp1, nksqr, lendiag, lenud, lenfull, nmax;
    mwIndex k, i, ii;
    double *lab, *q, *labk;
    const double *x;
    
    /* Argument check */
    mxAssert(nrhs >= NPARIN, "psdeig requires more input arguments");
    mxAssert(nlhs <= NPAROUT, "psdeig produces less output arguments");
    
    /* Disassemble cone structure and determine output sizes */
    conepars(K_IN, &cK);
    lendiag = cK.rLen + cK.hLen;
    lenud = cK.rDim + cK.hDim;
    lenfull = cK.lpN + cK.qDim + lenud;
    
    /* Get input vector x and skip to the PSD terms */
    mxAssert( !mxIsSparse(X_IN), "x must be full (not sparse)." );
    x = mxGetPr(X_IN);
    if ( mxGetM(X_IN) * mxGetN(X_IN) != lenud ) {
        mxAssert( mxGetM(X_IN) * mxGetN(X_IN) == lenfull, "Size mismatch x" );
        x += cK.lpN + cK.qDim;       
    }
    
    /* Allocate the output arrays */
    LAB_OUT = mxCreateDoubleMatrix( lendiag, (mwSize)1, mxREAL );
    lab = mxGetPr( LAB_OUT );
    if ( nlhs > 1 ) {
        Q_OUT = mxCreateDoubleMatrix( lenud, (mwSize)1, mxREAL );
        q = mxGetPr( Q_OUT );
    }
  
    /* 
     * Real symmetric matrices   
     */
    if ( cK.rsdpN != 0 ) {
        nmax = 1;
        for ( k = 0 ; k != cK.rsdpN ; ++k ) {
            nk = (mwSize)cK.sdpNL[k];
            if ( nmax < nk ) nmax = nk;
        }
        Xk = mxCreateDoubleMatrix( nmax, nmax, mxREAL );
        for ( k = 0 ; k != cK.rsdpN ; ++k ) {
            nk = (mwSize)cK.sdpNL[k]; 
            nksqr = SQR(nk);
            mxSetM( Xk, nk ); mxSetN( Xk, nk );
            /* Symmetric projection onto the work array */
            symproj( mxGetPr(Xk), x, nk );
            if ( nlhs <= 1 ) {
                /* One argument only: lab */
                mexCallMATLAB( 1, output_array, 1, &Xk, "eig" );
                memcpy( lab, mxGetPr(output_array[0]), nk * sizeof(double) );
                mxDestroyArray( output_array[0] );
            } else {
                /* First argument: Q */
                mexCallMATLAB( 2, output_array, 1, &Xk, "eig" );
                memcpy( q, mxGetPr(output_array[0]), nksqr * sizeof(double) );
                q += nksqr;
                /* Second argument: extract diag(Lab) */
                nkp1 = nk + 1;
                labk = mxGetPr( output_array[1] );
                for(i = 0, ii = 0; i < nk; i++, ii += nkp1)
                    lab[i] = labk[ii];
                mxDestroyArray( output_array[0] );
                mxDestroyArray( output_array[1] );
            }
            lab += nk;
            x += nksqr;
        }
        mxDestroyArray( Xk );
    }
    
    /* Complex Hermitian matrices 
     * Note that if a complex argument is offered, even the "real" matrices
     * must be passed through here, just to be safe.
     */
    if ( cK.sdpN != cK.rsdpN ) {
        nmax = 1;
        for ( k = cK.rsdpN ; k != cK.sdpN ; ++k ) {
            nk = (mwSize)cK.sdpNL[k];
            if ( nmax < nk ) nmax = nk;
        }
        Xk = mxCreateDoubleMatrix( nmax, nmax, mxCOMPLEX );
        for ( k = cK.rsdpN ; k != cK.sdpN ; ++k ) {
            nk = (mwSize)cK.sdpNL[k]; 
            nksqr = SQR(nk);
            mxSetM(Xk, nk); mxSetN(Xk, nk);
            /* Skew-symmetric projection onto the work matrix */
            symproj( mxGetPr(Xk), x, nk );
            skewproj( mxGetPi(Xk), x + nksqr, nk );
            if ( nlhs <= 1 ) {
                /* One argument only: lab */
                mexCallMATLAB( 1, output_array, 1, &Xk, "eig" );
                memcpy(lab, mxGetPr(output_array[0]), nk * sizeof(double));
                mxDestroyArray(output_array[0]);
            } else {
                #ifdef USE_SVD
                mexCallMATLAB(3, output_array, 1, &Xk, "svd");
                #else
                mexCallMATLAB(2, output_array, 1, &Xk, "eig");
                #endif
                /* First argument: Q */
                memcpy( q, mxGetPr(output_array[0]), nksqr * sizeof(double) );
                q += nksqr;
                if( mxIsComplex(output_array[0]) )     /* if any imaginary part */
                    memcpy( q, mxGetPi(output_array[0]), nksqr * sizeof(double) );
                q += nksqr;
                /* Second argument: extract diag(Lab) */
                nkp1 = nk + 1;
                labk = mxGetPr(output_array[1]);
                for(i = 0, ii = 0; i < nk; i++, ii += nkp1)
                    lab[i] = labk[ii];
                mxDestroyArray(output_array[0]);
                mxDestroyArray(output_array[1]);
                #ifdef USE_SVD
                mxDestroyArray(output_array[2]);
                #endif
            }
            lab += nk;
            x += 2 * nksqr;
        }
        mxDestroyArray( Xk );
    }
    
}
