#include "mex.h"
#include <stddef.h>
#include <math.h>

#if !defined(HAVE_OCTAVE) && ( !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 ) )
typedef int mwIndex;
typedef int mwSize;
#endif

/*
% Classifications:
% 0  - zero
% 1  - negative constant
% 2  - real constant
% 3  - positive constant
% 4  - complex constant
% 5  - negative concave
% 6  - concave
% 7  - positive concave
% 8  - negative affine
% 9  - real affine
% 10 - positive affine
% 11 - negative convex
% 12 - convex
% 13 - positive convex
% 14 - complex affine
% 15 - log concave
% 16 - log affine
% 17 - log convex monomial
% 18 - log convex posynomial
% 19 - invalid
*/

static const char negmap[20] =
    { 0, 3, 2, 1, 4, 7, 6, 5,10, 9, 8,13,12,11,14,19,12, 6, 6,19 }; /* -x */
    
static const char cplxmap[20] =
    { 0, 4, 4, 4, 4,19,19,19,14,14,14,19,19,19,19,19,19,19,19,19 }; /* complex * x */
    
static const char addmap[20][20] = {
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19 }, /* zero + x */
    { 1, 1, 2, 2, 4, 5, 6, 6, 8, 9, 9,11,12,12,14,19,12,12,12,19 }, /* negative + x */
    { 2, 2, 2, 2, 4, 6, 6, 6, 9, 9, 9,12,12,12,14,19,12,12,12,19 }, /* real + x */
    { 3, 2, 2, 3, 4, 6, 6, 7, 9, 9,10,12,12,13,14,19,18,18,18,19 }, /* positive + x */
    { 4, 4, 4, 4, 4,19,19,19,14,14,14,19,19,19,14,19,19,19,19,19 }, /* complex + x */
    { 5, 5, 6, 6,19, 5, 6, 6, 5, 6, 6,19,19,19,19,19,19,19,19,19 }, /* n_concave + x */
    { 6, 6, 6, 6,19, 6, 6, 6, 6, 6, 6,19,19,19,19,19,19,19,19,19 }, /* concave + x */
    { 7, 6, 6, 7,19, 6, 6, 7, 6, 6, 7,19,19,19,19,19,19,19,19,19 }, /* p_concave + x */
    { 8, 8, 9, 9,14, 5, 6, 6, 8, 9, 9,11,12,12,14,19,12,12,12,19 }, /* n_affine + x */
    { 9, 9, 9, 9,14, 6, 6, 6, 9, 9, 9,12,12,12,14,19,12,12,12,19 }, /* affine + x */
    {10, 9, 9,10,14, 6, 6, 7, 9, 9,10,12,12,13,14,19,13,13,13,19 }, /* p_affine + x */
    {11,11,12,12,19,19,19,19,11,12,12,11,12,12,14,19,12,12,12,19 }, /* n_convex + x */
    {12,12,12,12,19,19,19,19,12,12,12,12,12,12,14,19,12,12,12,19 }, /* convex + x */
    {13,12,12,13,19,19,19,19,12,12,13,12,12,13,14,19,13,13,13,19 }, /* p_convex + x */
    {14,14,14,14,14,19,19,19,14,14,14,19,19,19,14,19,19,19,19,19 }, /* c_affine + x */
    {15,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19 }, /* l_concave + x */
    {16,12,12,18,19,19,19,19,12,12,13,12,12,13,19,19,18,18,18,19 }, /* l_affine + x */
    {17,12,12,18,19,19,19,19,12,12,13,12,12,13,19,19,18,18,18,19 }, /* monomial + x */
    {18,12,12,18,19,19,19,19,12,12,13,12,12,13,19,19,18,18,18,19 }, /* posynomial + x */
    {19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19 } };

/*
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
*/

void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )

{
    if ( !mxIsSparse( prhs[0] ) )
        mexErrMsgTxt( "Requires a sparse argument" );
    mwIndex n = mxGetN( prhs[0] );
    plhs[0] = mxCreateNumericMatrix( (mwSize)1, n, mxINT8_CLASS, mxREAL );
    if ( plhs[0] == 0 )
        mexErrMsgTxt( "Unable to allocate output argument" );
    mwIndex  mB   = mxGetNumberOfElements( prhs[1] );
    const mwIndex* ir   = mxGetIr( prhs[0] );
    const mwIndex* jc   = mxGetJc( prhs[0] );
    const double*  pr   = mxGetPr( prhs[0] );
    const double*  pi   = mxGetPi( prhs[0] );
    const char*    vexs = (char*)mxGetData( prhs[1] );
    char*    ans  = (char*)mxGetData( plhs[0] );
    double   ival = 0;
    for ( mwIndex c = 0 ; c != n ; ++c, ++jc ) {
        int vex = 0;
        for ( mwIndex m = jc[0], mEnd = jc[1] ; m != mEnd ; ++m ) {
            double p = pr[m];
            if ( p == 0 ) continue;
            mwIndex r = ir[m];
            int nvex = r >= mB ? 9 : vexs[r];
            if ( isnan(p) )
                nvex = 19;
            else if ( isinf(p) ) {
                if ( nvex == 3 ) ival += p;
                else if ( nvex == 1 ) ival -= p;
                else nvex = 19;
            } else if ( p < 0 )
                nvex = negmap[nvex];
            vex = addmap[vex][nvex];
            if ( vex == 19 ) break;
        }
        if ( pi != 0 && vex != 19 ) {
            for ( mwIndex m = jc[0], mEnd = jc[1] ; m != mEnd ; ++m ) {
                double p = pi[m];
                if ( p == 0 ) continue;
                if ( isnan(p) || isinf(p) ) { vex = 19; break; }
                mwIndex r = ir[m];
                int nvex = r >= mB ? 14 : cplxmap[(int)vexs[r]];
                vex = addmap[vex][nvex];
                if ( vex == 19 ) break;
            }
        }
        if ( ival != 0 && vex != 19 ) {
            if ( isnan(ival) ) vex = 19;
            else vex = ival < 0 ? 1 : 3;
            ival = 0;
        }
        ans[c] = vex ? vex : 2;
    }
}
