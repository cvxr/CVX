function [ x, y ] = cvx_sparse_if( x, y )

if nargin < 2,
    y = x;
end
if islogical( y ),
    tst = y;
else
    [ m, n ] = size( y );
    d = 2 * ( 1 + ~isreal( y ) );
    tst = ( d + 1 ) * nnz( y ) < ( d * m - 1 ) * n;
end
if tst
    x = sparse( x );
    if nargout > 2,
        y = sparse( y );
    end
else
    x = full( x );
    if nargout > 2,
        y = full( y );
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
    
