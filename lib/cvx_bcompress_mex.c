#include "mex.h"
#include "matrix.h"

/*
% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
*/

#define MY_EPS 1.776356839400251e-015
int *ss_ir, *ss_jc, *ss_temp, *ss_ndxs, ss_sign;
double *ss_pr;

static int mycomp( int c1, int c2 )
{
    int d1, d2;
    int *j1, *j2, *rb1, *rb2, *re1, *re2;
    double nm1, nm2, b1, b2, bs, bd;
    double *pb1, *pb2, *pe1;

    /*
     * Comparison #1: number of elements
     */

    j1 = ss_jc + c1;
    j2 = ss_jc + c2;
    d1 = j1[1] - j1[0];
    d2 = j2[1] - j2[0];
    if ( d1 < d2 ) return -1;
    if ( d1 > d2 ) return +1;
    if ( d1 == 0 ) return  0;
    
    /*
     * Comparison #2: Lexical pattern
     */

    rb1 = ss_ir + j1[0];
    re1 = rb1   + d1;
    rb2 = ss_ir + j2[0];
    re2 = rb2   + d2;
    while ( 1 ) {
        if ( *rb1 < *rb2  ) return -1;
        if ( *rb1 > *rb2  ) return +1;
        if ( ++rb1 == re1 )
            if ( ++rb2 == re2 ) break;
            else return -1;
        else if ( ++rb2 == re2 )
            return +1;
    }
    
    /*
     * Comparison #3: sort by normalized column values
     */

    pb1 = ss_pr + j1[0];
    pb2 = ss_pr + j2[0];
    /* Quick exit for one-element columns */
    if ( d1 == 1 )
        if ( ss_sign )       return 0;
        else if ( *pb1 < 0 ) return *pb2 < 0 ?  0 : -1;
        else                 return *pb2 < 0 ? +1 :  0;
    pe1 = pb1 + d1;
    nm1 = 1.0 / *pb1++;
    nm2 = 1.0 / *pb2++;
    if ( !ss_sign ) {
        if ( nm1 < 0 ) nm1 = -nm1;
        if ( nm2 < 0 ) nm2 = -nm2;
    }
    while ( 1 ) {
        b1 = nm1 * *pb1;
        b2 = nm2 * *pb2;
        bd = b1 - b2;
        bs = b1 + b2;
        if ( ( bd < 0 ? -bd : +bd ) > MY_EPS * ( bs < 0 ? -bs : +bs ) )
            return b1 < b2 ? -1 : +1;
        if ( ++pb1 == pe1 ) break;
        ++pb2;
    }
   
    return 0;
}

void merge( int *nb, int* nm, int* ne )
{
    int *i = nb, *j = nm, *t = ss_temp;
    for ( ;; )
        if ( mycomp( *i, *j ) <= 0 ) {
            *t++ = *i++;
            if ( i == nm ) {
                while ( j != ne )
                    *t++ = *j++;
                break;
            }
        } else {
            *t++ = *j++;
            if ( j == ne ) {
                while ( i != nm )
                    *t++ = *i++;
                break;
            }
        }
    for ( i = nb, t = ss_temp ; i != ne ; *i++ = *t++ );
}

void merge_sort( int *nb, int* ne )
{
    int *nm;
    ptrdiff_t numel = ne - nb;
    if ( numel < 2 ) return;
    nm = nb + ( numel / 2 );
    merge_sort( nb, nm );
    merge_sort( nm, ne );
    merge( nb, nm, ne );
}

void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )

{
    static const char* errmsg = "Error allocating temporary data";
    int n, n2, k, col, col2, nsrt, *ss_ndxs, *p, lastcol, usesign;
    double norm, sign, *map, *scl, lastnorm;
    
    n       = mxGetN(  prhs[0] );
    ss_ir   = mxGetIr( prhs[0] );
    ss_jc   = mxGetJc( prhs[0] );
    ss_pr   = mxGetPr( prhs[0] );
    ss_sign = nrhs < 2 || mxGetScalar(prhs[1]) == 0;
    nsrt    = nrhs < 3 ? 0 : (int)mxGetScalar(prhs[2]);
    ss_ndxs = mxCalloc( n, sizeof(int) );
    ss_temp = mxCalloc( n, sizeof(int) );
    
    /*
     * Sort the column indices
     */

    for ( col = 0 ; col != n ; ++col )
        ss_ndxs[col] = col;
    if ( nsrt < n ) {
        merge_sort( ss_ndxs + nsrt, ss_ndxs + n );
        if ( nsrt > 0 )
            merge( ss_ndxs, ss_ndxs + nsrt, ss_ndxs + n );
    }

    /*
     * Determine which rows are unique, and which are scales of another row
     */

    plhs[0] = mxCreateDoubleMatrix( 1, n, mxREAL );
    plhs[1] = mxCreateDoubleMatrix( 1, n, mxREAL );
    if ( plhs[0] == 0 || plhs[1] == 0 )
        mexErrMsgTxt( "Unable to allocate output arguments" );
    map = mxGetPr( plhs[0] );
    scl = mxGetPr( plhs[1] );
    lastcol = -1;
    lastnorm = 0.0;
    for ( k = 0 ; k != n ; ++k ) {
        col = ss_ndxs[k];
        if ( ss_jc[col] == ss_jc[col+1] ) {
            map[col] = col + 1;
            scl[col] = 0;
        } else if ( lastcol == -1 || mycomp( lastcol, col ) != 0 ) {
            lastcol  = col;
            lastnorm = ss_pr[ss_jc[col]];
            map[col] = lastcol + 1;
            scl[col] = 1.0;
        } else {
            map[col] = lastcol + 1;
            scl[col] = lastnorm / ss_pr[ss_jc[col]];
        }
    }
}
