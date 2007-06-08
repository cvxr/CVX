function y = cvx_recip( x )

%CVX_RECIP    Reciprocal with a check for divide by zero.
%
%   CVX_RECIP(X) is the elementwise inverse of X; i.e., 1.0 ./ X.
%       Unlike the native MATLAB division, division by zero is NOT
%       permitted and will result in an error.
%
%   Disciplined quadratic programming information:
%       When used in CVX expressions, X must be constant, monomial,
%       or posynomial.

error( nargchk( 1, 1, nargin ) );
if any( x( : ) == 0 ),
    error( sprintf( 'Disciplined convex programming error:\n    Division by zero.' ) );
else
    y = 1.0 ./ x;
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
