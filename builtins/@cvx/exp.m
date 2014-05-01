function y = exp( x )

%   Disciplined convex programming information:
%       EXP(X) is convex and nondecreasing in X. When used in CVX
%       expressions, X must be real. Typically, X must also be affine
%       or convex; X can also be concave, but this produces a log-concave
%       result with very limited usefulness.
%
%   Disciplined geometric programming information:
%       EXP(X) is typically not used in geometric programs. However,
%       EXP(X), where X is a monomial or posynomial, can be included in 
%       geometric programs wherever a posynomial would be appropriate.

cvx_expert_check( 'exp', x );
try
    y = cvx( x.size_, cvx_getexp( x.basis_, true ) );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ) throw( exc );
    else rethrow( exc ); end
end
        
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
