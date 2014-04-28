function zD = cvx_pushcnstr( zb, iseq )

global cvx___
mO = cvx___.n_equality;
zi = ~isreal( zb );
if zi, zb = [ real(zb), imag(zb) ]; end
mN = size( zb, 2 );
cvx___.equalities{end+1} = zb;
cvx___.needslack(end+1,1) = ~iseq;
cvx___.n_equality = cvx___.n_equality + mN;
if nargout,
    zD = sparse( mO + 1 : mO + mN, 1 : mN, 1 );
    if zi, zD = zD(:,1:end/2) + 1j * zD(:,end/2+1:end); end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
