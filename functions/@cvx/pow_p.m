function y = pow_p( x, p )

%POW_P   Internal cvx version.

error(nargchk(2,2,nargin)); %#ok
y = pow_cvx( x, p, 'pow_p' );

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
