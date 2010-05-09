function y = cvx_isnonzero( x, full )
error( nargchk( 1, 2, nargin ) );
y = any( x.basis_, 1 );
if nargin < 2,
    y = all( y );
else
    y = cvx_reshape( y, x.size_ );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
