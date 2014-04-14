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
    mneg = int8([1,2,2,19,4,4,2,8,8,2,11,11,2,19,19,19,19,19,19]);
end
cvx___.classes(tx) = mneg(cvx___.classes(tx));

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
