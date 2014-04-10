#include "mex.h"
#include <stddef.h>
#include <math.h>

#if !defined(HAVE_OCTAVE) && ( !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 ) )
typedef int mwIndex;
typedef int mwSize;
#endif

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
    mwIndex  mA = mxGetM( prhs[0] );
    mwIndex  mB = mxGetNumberOfElements( prhs[1] );
    mwIndex  n  = mxGetN(  prhs[0] );
    mwIndex* ir = mxGetIr( prhs[0] );
    mwIndex* jc = mxGetJc( prhs[0] );
    double*  pr = mxGetPr( prhs[0] );
    double*  pi = mxGetPi( prhs[0] );
    double*  vexs = mxGetPr( prhs[1] );
    plhs[0] = mxCreateDoubleMatrix( (mwSize)1, n, mxREAL );
    if ( plhs[0] == 0 )
        mexErrMsgTxt( "Unable to allocate output arguments" );
    double* ans = mxGetPr( plhs[0] );
    
    if ( !mxIsSparse( prhs[0] ) ) {
        for ( mwIndex c = 0 ; c != n ; ++c, pr += mA ) {
            int vex = 0;
            for ( mwIndex r = 0 ; r != mA ; ++r ) {
                double p = pr[r];
                if ( p == 0 ) continue;
                if ( isnan(p) || ( r && isinf(p) ) ) { vex = 2; break; }
                if ( r >= mB ) break;
                int nvex = vexs[r];
                if ( nvex == 0 ) continue;
                if ( p < 0 ) nvex = -nvex;
                if ( vex == 0 ) vex = nvex;
                else if ( vex != nvex ) { vex = 2; break; }
            }
            ans[c] = vex == 2 ? NAN : vex;
        }
        if ( pi != 0 ) {
            for ( mwIndex c = 0 ; c != n ; ++c, pi += mA ) {
                int vex = isnan(ans[c]) ? 2 : (int)ans[c];
                if ( vex == 2 ) continue;
                for ( mwIndex r = 0 ; r != mA ; ++r ) {
                    double p = pi[r];
                    if ( p == 0 ) continue;
                    if ( vex || isnan(p) || ( r && isinf(p) ) ) { vex = 2; break; }
                    if ( r >= mB ) break;
                    if ( vexs[r] != 0 ) { vex = 2; break; }
                }
                if ( vex == 2 ) ans[c] = NAN;
            }
        }
    } else {
        for ( mwIndex c = 0, m = 0 ; c != n ; ++c, ++jc ) {
            int vex = 0;
            for ( mwIndex mEnd = jc[1] ; m != mEnd ; ++m ) {
                double p = pr[m];
                if ( p == 0 ) continue;
                mwIndex r = ir[m];
                if ( isnan(p) || ( r && isinf(p) ) ) { vex = 2; break; }
                if ( r >= mB ) break;
                int nvex = vexs[r];
                if ( nvex == 0 ) continue;
                if ( p < 0 ) nvex = -nvex;
                if ( vex == 0 ) vex = nvex;
                else if ( vex != nvex ) { vex = 2; break; }
            }
            ans[c] = vex == 2 ? nan(0) : vex;
        }
        if ( pi != 0 ) { 
            jc -= n;
            for ( mwIndex c = 0, m = 0 ; c != n ; ++c, ++jc ) {
                int vex = isnan(ans[c]) ? 2 : (int)ans[c];
                if ( vex == 2 ) continue;
                for ( mwIndex mEnd = jc[1] ; m != mEnd ; ++m ) {
                    double p = pi[m];
                    if ( pi == 0 ) continue;
                    if ( vex || isnan(p) ) { vex = 2; break; }
                    mwIndex r = ir[m];
                    if ( r && isinf(p) ) { vex = 2; break; }
                    if ( r >= mB ) break;
                    if ( vexs[r] != 0 ) { vex = 2; break; }
                }
            }
        }
    }
    
}
