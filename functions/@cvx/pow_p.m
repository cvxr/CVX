function y = pow_p( x, p )

%POW_P   Internal cvx version.

error(nargchk(2,2,nargin));
y = pow_cvx( x, p, 'pow_p' );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
