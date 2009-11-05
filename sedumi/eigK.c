/*
   [lab,q] = eigK(x,K)
   Computes spectral coefficients of x w.r.t. K
   Arguments "q" is optional - without it's considerably
   faster in case of PSD blocks.
   FLOPS indication: 1.3 nk^3 versus 9.0 nk^3 for nk=500,
                     1.5 nk^3        9.8 nk^3 for nk=50.
 
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

#define LAB_OUT plhs[0]
#define Q_OUT plhs[1]
#define NPAROUT 2

#define X_IN prhs[0]
#define K_IN prhs[1]
#define NPARIN 2

/*
 * QEIG: computes the 2 spectral values w.r.t. Lorentz cone
 * INPUT: x, n - vector of length n > 1
 * OUTPUT: lab - length-2 vector of spectral values
 */

void qeig( double *lab, const double *x, const mwIndex n )
{
    double nx2;
    nx2    = sqrt( realssqr( x+1, n-1 ) );
    lab[0] = ( x[0] - nx2 ) / M_SQRT2;
    lab[1] = ( x[0] + nx2 ) / M_SQRT2;
}

/*
 * CXQEIG: computes the 2 spectral values w.r.t. complex Lorentz cone
 * INPUT: x, xpi, n - real and imaginary parts of length n > 1
 * OUTPUT: lab - length-2 vector of spectral values
 */

void cxqeig( double *lab, const double *x, const double *xpi, const mwSize n )
{
    double nx2;
    nx2    = sqrt( realssqr( x+1, n-1 ) + realssqr( xpi+1, n-1 ) );
    lab[0] = ( x[0] - nx2 ) / M_SQRT2;
    lab[1] = ( x[0] + nx2 ) / M_SQRT2;
}

/*
 * RCONEIG: computes the 2 spectral values w.r.t. rotated Lorentz cone
 * INPUT: x, n - vector of length n > 2
 * OUTPUT: lab - length-2 vector of spectral values
 */

void rconeeig_i( double *lab, double x1, double x2, double x3sqr )
{
    double t1, t2, t3;
    t1 = x1 + x2;
    t2 = x1 * x2 - x3sqr / 2;
    t3 = sqrt( SQR(t1) - 4 * t2 );
    lab[0] = ( t1 + ( t1 < 0 ? -t3 : t3 ) ) / 2;
    lab[1] = lab[0] == 0 ? t1 : t2 / lab[0];
}

void rconeeig( double *lab, const double* x, const mwSize n )
{
    rconeeig_i( lab, x[0], x[1], realssqr(x+2,n-2) );
}

/*
 * CXRCONEIG: computes the 2 spectral values w.r.t. rotated Lorentz cone
 * INPUT: xpi, n - complex vector of length n > 2
 * OUTPUT: lab - length-2 vector of spectral values
 */

void cxrconeeig( double *lab, const double* x, const double* xpi, const mwSize n )
{
    /* roots of: lab^2 - (x1+x2)*lab + (x1*x2-x3sqr/2) = 0 */
    rconeeig_i( lab, x[0], x[1], realssqr(x+2,n-2) + realssqr(xpi+2,n-2) );
}

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
 *   [lab,q] = eigK(x,K)
 *   Computes spectral coefficients of x w.r.t. K
 * REMARK
 *    If this function is used internally by SeDuMi, then
 *    complex numbers are stored in a single real vector. To make
 *    it invokable from the Matlab command-line by the user, we
 *    also allow Matlab complex vector x.
 */

void mexFunction(
const int nlhs, mxArray *plhs[],
const int nrhs, const mxArray *prhs[])
{
    mxArray *output_array[3], *Xk;
    coneK cK;
    int xcplx;
    mwSize nk, nksqr, lendiag, lenfull, nmax, nkp1;
    mwIndex k, i, ii;
    double *lab, *q, *qpi, *labk;
    const double *x, *xpi;
    
    /* Argument check */
    mxAssert(nrhs >= NPARIN, "psdeig requires more input arguments");
    mxAssert(nlhs <= NPAROUT, "psdeig produces less output arguments");
    
    /* Disassemble cone structure and determine output sizes */
    conepars(K_IN, &cK);
    lendiag = cK.lpN + 2 * (cK.lorN + cK.rconeN) + cK.rLen + cK.hLen;
    lenfull = cK.lpN + cK.qDim + cK.rDim + cK.hDim;
    if( cK.rconeN > 0 )
        for( i = 0 ; i != cK.rconeN ; ++i )
            lenfull += (mwSize)cK.rconeNL[i];

    /* Get input vector x */
    mxAssert( mxGetM(X_IN) * mxGetN(X_IN) == lenfull, "Size mismatch x" );
    mxAssert( !mxIsSparse(X_IN), "x must be full (not sparse)." );
    x = mxGetPr(X_IN);
    xcplx = mxIsComplex(X_IN);
    if ( xcplx )
        xpi = mxGetPi(X_IN);
    
    /* Allocate the output arrays */
    LAB_OUT = mxCreateDoubleMatrix( lendiag, (mwSize)1, mxREAL );
    lab = mxGetPr( LAB_OUT );
    if ( nlhs > 1 ) {
        if ( mxIsComplex(X_IN) ) 
            Q_OUT = mxCreateDoubleMatrix( cK.rDim, (mwSize)1, mxCOMPLEX );
        else
            Q_OUT = mxCreateDoubleMatrix( cK.rDim + cK.hDim, (mwSize)1, mxREAL );
        q = mxGetPr( Q_OUT );
        qpi = mxGetPi( Q_OUT );
    }
    
    /* Linear terms: lab = x */
    if ( cK.lpN ) {
        memcpy( lab, x, cK.lpN * sizeof(double) );
        lab += cK.lpN;
        x += cK.lpN;
        if ( xcplx )
            xpi += cK.lpN;
    }
    
    /* Lorentz cones  */
    for ( k = 0 ; k != cK.lorN ; ++k ) {
        nk = (mwSize)cK.lorNL[k];
        if ( xcplx ) {
            cxqeig( lab, x, xpi, nk );
            xpi += nk;
        } else
            qeig( lab, x, nk );
        lab += 2;
        x += nk;
    }
    
    /* Rotated Lorentz cones */
    for ( k = 0 ; k != cK.rconeN ; ++k ) {
        nk = (mwSize)cK.rconeNL[k];
        if ( xcplx ) {
            cxrconeeig( lab, x, xpi, nk );
            xpi += nk;
        } else
            rconeeig( lab, x, nk );
        lab += 2;
        x += nk;
    }
    
    /* 
     * Real symmetric matrices   
     * Note that if a complex argument is offered, then all of the matrices
     * are treated as complex, so we have to skip this stage
     */
    if ( !xcplx && cK.rsdpN != 0 ) {
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
    ii = xcplx ? 0 : cK.rsdpN;
    if ( cK.sdpN != ii ) {
        nmax = 1;
        for ( k = ii ; k != cK.sdpN ; ++k ) {
            nk = cK.sdpNL[k];
            if ( nmax < nk ) nmax = nk;
        }
        Xk = mxCreateDoubleMatrix( nmax, nmax, mxCOMPLEX );
        for ( k = ii ; k != cK.sdpN ; ++k ) {
            nk = cK.sdpNL[k]; nksqr = SQR(nk);
            mxSetM(Xk, nk); mxSetN(Xk, nk);
            /* Skew-symmetric projection onto the work matrix */
            symproj( mxGetPr(Xk), x, nk );
            skewproj( mxGetPi(Xk), xcplx ? xpi : x + nksqr, nk );
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
            if ( xcplx ) {
                x += nksqr;
                xpi += nksqr;
            } else
                x += 2 * nksqr;
        }
        mxDestroyArray( Xk );
    }
    
}
