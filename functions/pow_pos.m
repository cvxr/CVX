function y = pow_pos( x, p )

%POW_POS   Power of positive part.
%   POW_POS(X,P) = MAX(X,0).^P. P must be real.
%
%   Disciplined convex programming information:
%       The geometry of POW_POS(X,P) depends on the precise value of P,
%       which must be a real constant:
%                P < 0: convex and nonincreasing;  X must be concave.
%           0 <= P < 1: concave and nondecreasing; X must be concave.
%           1 <= P    : convex and nondecreasing;  X must be convex.
%       In all cases, X must be real.

error( nargchk( 2, 2, nargin ) );
if ~isnumeric( p ) | ~isreal( p ),
    error( 'P must be real.' );
elseif ~isreal( x ),
    error( 'X must be real.' );
end
y = max(x,0).^p;

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
