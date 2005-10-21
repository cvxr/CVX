#include "mex.h"
#include "matrix.h"
#include "assert.h"

/*
% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
*/

typedef struct {
    int *values;
    int *temp;
} col_sort_struct;

static int mycomp( int i, int j, col_sort_struct* ss )
{
    int vi = ss->values[i], vj = ss->values[j];
    return vi > vj ? -1 : ( vi < vj ? +1 : ( i < j ? -1 : +1 ) );
}

static void merge( int *nb, int* nm, int* ne, col_sort_struct* ss )
{
    int *i = nb, *j = nm, *t = ss->temp;
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

void merge_sort( int *nb, int* ne, col_sort_struct* ss )
{
    int *nm;
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
    int m            =  mxGetM( prhs[0] ),
        n            =  mxGetN( prhs[0] ),
        *row_index_A = mxGetIr( prhs[0] ),
        *col_count_A = mxGetJc( prhs[0] ),
        nnzA         = col_count_A[n],     
        nnzA1        = nnzA + 1,
        *row_index_C = mxGetIr( prhs[1] ),
        *col_count_C = mxGetJc( prhs[1] ),
        *row_counts  = (int*)mxCalloc( m + n, sizeof(int) ),
        *candidates  = row_counts + m,
        *row_flags   = (int*)mxGetData( plhs[0] = mxCreateNumericMatrix( m, 1, mxINT32_CLASS, mxREAL ) ),
        *col_flags   = (int*)mxGetData( plhs[1] = mxCreateNumericMatrix( 1, n, mxINT32_CLASS, mxREAL ) ),
        *reserved    = (int*)mxGetData( prhs[2] ),
        *nreserved   = (int*)mxGetData( plhs[2] = mxCreateNumericMatrix( 1, n, mxINT32_CLASS, mxREAL ) ),
        j, nTry, nQ, nnzQ, pass;
    assert( sizeof(int) == 4 );
    
    for ( j = 0 ; j != nnzA ; ++j )
        ++row_counts[row_index_A[j]];
    nQ = nnzQ = 0;
    for ( j = 0 ; j != n ; ++j )
        if ( ( nreserved[j] = reserved[j] ) == 0 ) {
            int kBeg = col_count_A[j],
                kEnd = col_count_A[j+1], k,
                nc   = ( kEnd - kBeg ) + ( col_count_C[j+1] - col_count_C[j] ) - 1,
                found = 0, rr = nnzA1, rb;
            for ( k = kBeg ; k != kEnd ; ++k ) {
                int tr = row_index_A[k],
                    nr = row_counts[tr] - 1,
                    n1 = nr * nc,
                    n2 = nr + nc + 1;
                if ( n1 <= n2 ) {
                    found = 1;
                    if ( row_flags[tr] <= 0 && ( n1 < rr || n1 == rr && row_flags[rb] && !row_flags[tr] ) ) {
                        rb = tr;
                        rr = n1;
                    }
                }
            }
            if ( rr != nnzA1 ) {
                for ( k = kBeg ; k != kEnd ; ++k ) {
                    int tr = row_index_A[k],
                        fg = row_flags[tr];
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
            } else if ( !found )
                nreserved[j] = 1;
        }
    for ( j = 0 ; j != m ; ++j )
        if ( row_flags[j] < 0 )
            row_flags[j] = 0;
    if ( nnzQ > nQ ) {
        col_sort_struct ss;
        int i, iEnd;
        ss.values = (int*)mxCalloc( 2 * n, sizeof(int) );
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
            } else
                nreserved[j] = 1;
        }
        nnzQ /= 2;
        merge_sort( candidates, candidates + iEnd, &ss );
        for ( i = 0 ; i != iEnd && nnzQ > 0 ; ++i ) {
            int j    = candidates[i],
                r    = col_flags[j] - 1,
                kBeg = col_count_A[j],
                kEnd = col_count_A[j+1], k;
            nnzQ -= row_flags[r] - 1;
            row_flags[r] = col_flags[j] = 0;
            for ( k = kBeg ; k != kEnd ; ++k )
                if ( row_flags[row_index_A[k]] != 0 )
                    --nnzQ;
        }
    }
    for ( j = 0 ; j != nQ ; ++j ) {
        int tc = candidates[j];
        if ( col_flags[tc] != 0 )
            nreserved[tc] = 1;
    }
}
