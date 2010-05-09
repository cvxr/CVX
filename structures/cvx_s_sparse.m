function y = cvx_s_sparse( m, n, i, j )
%CVX_S_SPARSE Matrices with a fixed sparsity pattern.

if nargin < 3,
    error( 'Sparsity structure missing.' );
elseif ~isnumeric( i ) || ~isnumeric( j ),
    error( 'Sparsity arguments must be vectors of nonnegative integers.' );
end
i = i(:);
j = j(:);
if any( i <= 0 ) || any( j <= 0 ) || any( i ~= floor( i ) ) || any( j ~= floor( j ) ),
    error( 'Sparsity arguments must be vectors nonnegative integers.' );
elseif length( i ) ~= 1 && length( j ) ~= 1 && length( i ) ~= length( j ),
    error( 'Sparsity arguments have incompatible size.' );
elseif any( i > m ) || any( j > n ),
    error( 'One or more indices are out of range.' );
end
nz = max( length( i ), length( j ) );
y = min( sparse( 1 : nz, i + ( j - 1 ) * m, 1, nz, m * n ), 1 );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.



