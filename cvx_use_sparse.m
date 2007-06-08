function tf = cvx_use_sparse( sz, nz, isr )

% CVX_USE_SPARSE    Sparse matrix size query
% 
% CVX_USE_SPARSE(SZ,NZ,ISR) is an internal function used by CVX to determine if
% a matrix of size SZ with NZ nonzeros would require less memory to store in
% sparse form or dense form. The matrix is real if ISR=1, and complex otherwise.

if nargin == 1,
    ss = size( sz );
    if any( ss == 1 ) | length( ss ) > 2,
        tf = false;
        return
    end
    isr = isreal( sz );
    if issparse( sz ),
        nz = nzmax( sz );
    else
        nz = nnz( sz );
    end
    sz = ss;
elseif any( sz == 1 ) | length( sz ) > 2,
    tf = false;
    return
elseif nargin < 3,
    isr = true;
end
if isr,
    tf = 1 + ( 1 - 2 * sz( 1 ) ) * sz( 2 ) + 3 * nz < 0;
else
    tf = 1 + ( 1 - 4 * sz( 1 ) ) * sz( 2 ) + 5 * nz < 0;
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
