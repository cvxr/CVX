function y = pow_abs( x, p )

%POW_POS   Power of absolute value.
%   POW_POS(X,P) = ABS(X).^P. P must be greater than or equal to one.
%
%   Disciplined convex programming information:
%       POW_ABS(X,P) is convex and nonmonotonic, so X must be affine.
%       P must be constant, and its elements must be greater than or
%       equal to one. X may be complex.

error( nargchk( 2, 2, nargin ) );
if ~isnumeric( p ) | ~isreal( p ),
    error( 'P must be real.' );
elseif any( p(:) < 1 ),
    error( 'P must be greater than or equal to one.' );
end
y = abs(x).^p;

% Copyright 2008 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
