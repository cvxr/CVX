function cvx_setnpos( x )

global cvx___
x  = x.basis_;
xc = x(1,:);
x  = x(2:end,:);
tt = ( xc <= 0 ) & ( sum( x ~= 0 ) == 1 );
[ tx, dummy ] = find( x(:,tt) ); %#ok
tx = tx + 1;
persistent mneg;
if isempty( mneg ),
    mneg = int8([1,2,2,22,5,5,6,8,8,9,11,11,12,22,22,22,22,22,22,22,22,22]);
end
cvx___.classes(tx) = mneg(cvx___.classes(tx));

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
