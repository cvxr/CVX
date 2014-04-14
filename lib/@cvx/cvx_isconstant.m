function y = cvx_isconstant( x, fz )
b = x.basis_;
if size( b, 1 ) <= 1,
    if nargin < 2 || ~fz,
    	y = true;
    else
    	y = true( size( x ) );
    end
elseif nargin < 2 || ~fz,
    y = full( any( b, 2 ) );
    y( 1 ) = false;
    y = ~any( y );
else
	y = any( b( 2 : end, : ), 1 );
    y = ~reshape( fz( y ), x.size_ );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
