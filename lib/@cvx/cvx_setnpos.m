function cvx_setnpos( x )

global cvx___
x  = x.basis_;
xc = x(1,:);
x  = x(2:end,:);
tt = ( xc <= 0 ) & ( sum( x ~= 0 ) == 1 );
[ tx, dummy ] = find( x(:,tt) ); %#ok
cvx___.sign(tx+1) = -1;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
