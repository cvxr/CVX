#include "mex.h"
#include "matrix.h"

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

/*
% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
*/

typedef struct {
    mwSize *values;
    mwSize *temp;
} col_sort_struct;

static int mycomp( mwIndex i, mwIndex j, col_sort_struct* ss )
{
    mwSize vi = ss->values[i], vj = ss->values[j];
    return vi > vj ? -1 : ( vi < vj ? +1 : ( i < j ? -1 : +1 ) );
}

static void merge( mwSize *nb, mwSize* nm, mwSize* ne, col_sort_struct* ss )
{
    mwSize *i = nb, *j = nm, *t = ss->temp;
    for ( ;; )
        if ( mycomp( *i, *j, ss ) <= 0 ) {
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
    for ( i = nb, t = ss->temp ; i != ne ; *i++ = *t++ );
}

void merge_sort( mwSize *nb, mwSize* ne, col_sort_struct* ss )
{
    mwSize *nm;
    ptrdiff_t numel = ne - nb;
    if ( numel < 2 ) return;
    nm = nb + ( numel / 2 );
    merge_sort( nb, nm, ss );
    merge_sort( nm, ne, ss );
    merge( nb, nm, ne, ss );
}

void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )

{
    mwIndex *row_index_A = mxGetIr( prhs[0] ),
            *col_count_A = mxGetJc( prhs[0] );
    mwSize  m            =  mxGetM( prhs[0] ),
            n            =  mxGetN( prhs[0] ),
            nobj         = (mwSize)*mxGetPr( prhs[1] ),
            nnzA         = (mwSize)col_count_A[n],
            nnzA1        = nnzA + 1,
            j, nQ, nnzQ,
            *row_counts  = (mwSize*)mxCalloc( m + n, sizeof(mwSize) ),
            *candidates  = row_counts + m;
    double  *row_flags   = mxGetPr( plhs[0] = mxCreateNumericMatrix( m, 1, mxDOUBLE_CLASS, mxREAL ) ),
            *col_flags   = mxGetPr( plhs[1] = mxCreateNumericMatrix( 1, n, mxDOUBLE_CLASS, mxREAL ) ),
            *reserved    = mxGetPr( prhs[2] ),
            *c_reserved  = mxGetPr( prhs[3] );

    for ( j = 0 ; j != nnzA ; ++j )
        ++row_counts[row_index_A[j]];
    nQ = nnzQ = 0;
    for ( j = nobj ; j != n ; ++j ) {
        if ( c_reserved[j] == 0 ) {
            mwIndex kBeg = col_count_A[j],
                    kEnd = col_count_A[j+1],
                    nc   = kEnd - kBeg, k, rb;
            mwSize  rr   = nnzA1;
            for ( k = kBeg ; k != kEnd ; ++k ) {
                mwIndex tr = row_index_A[k];
                if ( reserved[tr] == 0 ) {
                    mwSize nr = row_counts[tr],
                           n1 = ( nr - 1 ) * ( nc - 1 ),
                           n2 = nr + nc - 1;
                    if ( n1 <= n2 ) {
                        if ( row_flags[tr] <= 0 && ( n1 < rr || n1 == rr && row_flags[rb] && !row_flags[tr] ) ) {
                            rb = tr;
                            rr = n1;
                        }
                    }
                }
            }
            if ( rr != nnzA1 ) {
                for ( k = kBeg ; k != kEnd ; ++k ) {
                    mwIndex tr = row_index_A[k];
                    int32_T fg = row_flags[tr];
                    if ( fg > 0 ) {
                        ++fg;
                        ++nnzQ;
                    } else
                        --fg;
                    row_flags[tr] = fg;
                }
                col_flags[j] = rb + 1;
                nnzQ += ( row_flags[rb] = - row_flags[rb] );
                candidates[nQ++] = j;
            }
        }
    }
    for ( j = 0 ; j != m ; ++j )
        if ( row_flags[j] < 0 )
            row_flags[j] = 0;
    if ( nnzQ > nQ ) {
        col_sort_struct ss;
        int i, iEnd;
        ss.values = (mwSize*)mxCalloc( 2 * n, sizeof(mwSize) );
        ss.temp   = ss.values + n;
        nnzQ = 0;
        for ( i = iEnd = 0 ; i != nQ ; ++i ) {
            int j = candidates[i],
                r = col_flags[j] - 1,
                kBeg = col_count_A[j],
                kEnd = col_count_A[j+1], k;
            ss.values[j] = row_flags[r] - 2;
            for ( k = kBeg ; k != kEnd ; ++k )
                if ( row_flags[row_index_A[k]] != 0 )
                    ++ss.values[j];
            if ( ss.values[j] != 0 ) {
                candidates[iEnd++] = j;
                nnzQ += ss.values[j];
            }
        }
        nnzQ /= 2;
        merge_sort( candidates, candidates + iEnd, &ss );
        for ( i = 0 ; i != iEnd && nnzQ > 0 ; ++i ) {
            mwSize   j = candidates[i];
            mwIndex  r    = col_flags[j] - 1,
                     kBeg = col_count_A[j],
                     kEnd = col_count_A[j+1], k;
            mwSize  found = 0;
            for ( k = kBeg ; k != kEnd ; ++k ) {
                int tr = row_index_A[k];
                if ( r != tr && row_flags[tr] != 0 ) {
                    --row_flags[tr];
                    ++found;
                }
            }
            if ( found ) {
                nnzQ -= found + row_flags[r] - 1;
                row_flags[r] = col_flags[j] = 0;
            }
        }
    }
}
