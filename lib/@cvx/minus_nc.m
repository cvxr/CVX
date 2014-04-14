function z = minus_nc( x, y )

% Subtraction with no consistency checking.

x = cvx_basis( x );
y = cvx_basis( y );
mx = max(size(x),size(y));
x(end+1:mx(1),:) = 0;
y(end+1:mx(1),:) = 0;
z = bsxfun( @minus, x, y );
z = cvx( mx(2), z );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

