#include "mex.h"
#include <stddef.h>

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
    double*  sgns = mxGetPr( prhs[1] );
    double   dflt = mxGetPr( prhs[2] )[0];
    plhs[0] = mxCreateDoubleMatrix( (mwSize)1, n, mxREAL );
    if ( plhs[0] == 0 )
        mexErrMsgTxt( "Unable to allocate output arguments" );
    double* ans = mxGetPr( plhs[0] );
    
    if ( !mxIsSparse( prhs[0] ) ) {
        for ( mwIndex c = 0 ; c != n ; ++c, pr += mA, pi += mA ) {
            int first = 1, sign = dflt;
            for ( mwIndex r = 0 ; r != mA ; ++r ) {
                if ( pi && pi[r] ) { sign = 0; break; }
                if ( pr[r] == 0 ) continue;
                if ( r >= mB || sgns[r] == 0 ) { sign = 0; break; }
                int nsign = sgns[r] * pr[r] < 0 ? -1 : +1;
                if ( first ) { first = 0; sign = nsign; }
                else if ( sign != nsign ) { sign = 0; break; }
            }
            ans[c] = sign;
        }
    } else {
        for ( mwIndex c = 0, m = 0 ; c != n ; ++c, ++jc ) {
            int first = 1, sign = dflt;
            for ( mwIndex mEnd = jc[1] ; m != mEnd ; ++m ) {
                if ( pi && pi[m] ) { sign = 0; break; }
                if ( pr[m] == 0 ) continue;
                mwIndex r = ir[m];
                if ( r >= mB || sgns[r] == 0 ) { sign = 0; break; }
                int nsign = sgns[r] * pr[m] < 0 ? -1 : +1;
                if ( first ) { first = 0; sign = nsign; }
                else if ( sign != nsign ) { sign = 0; break; }
            }
            ans[c] = sign;
        }
    }
    
}
