function y = cvx_isnonzero( x, full )
error( nargchk( 1, 2, nargin ) );
y = any( x.basis_, 1 );
if nargin < 2,
    y = all( y );
else
    y = cvx_reshape( y, x.size_ );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
