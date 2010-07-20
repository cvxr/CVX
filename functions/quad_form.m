function y = quad_form( x, Q, v, w )

%QUAD_FORM quadratic form.
%   QUAD_FORM(x,Q) is x'*((Q+Q')/2)*x = real(x'*Q*x).
%   QUAD_FORM(x,Q,v,w) is real(x'*Q*x) + v(:)'*x + w.
%
%   x must be a row or column vector, and Q must either be a scalar or
%   a square matrix with the same number of rows as x(:).
%  
%   NOTE: The use of QUAD_FORM can often be replaced by a call to NORM. For
%   example, if Q is positive definite, then the constraint
%       quad_form(x,Q) <= 1
%   is equivalent to
%       norm(sqrtm(Q)*x) <= 1
%   Generally speaking, the NORM version will be more reliable and more
%   accurate, so we encourage you to make similar conversions whenever
%   possible. We *strongly* discourage the QP-era practice of converting
%   NORM expressions into quadratic forms. 
%
%   Disciplined convex programming information:
%       QUAD_FORM(x,Q) is neither convex nor concave in x and Q jointly,
%       so at least one of the two arguments must be constant.
%
%       If Q is constant, then QUAD_FORM is convex if Q is positive
%       semidefinite, and concave if Q is negative semidefinite. An error 
%       is generated if Q is indefinite (unless x is also constant). 
%       QUAD_FORM is nonmonotonic in x, so x must be affine.
%       
%       If x is constant, then QUAD_FORM is affine in Q, v, and w. The
%       signs of x will govern whether the elements of Q, v, and w may
%       be convex, concave, or affine.

error( nargchk( 2, 4, nargin ) );
sx = size( x );
if length( sx ) ~= 2 || all( sx > 1 ),
    error( 'The first argument must be a row or column.' );
elseif ndims( Q ) > 2 || size( Q, 1 ) ~= size( Q, 2 ),
    error( 'The second argument must be a scalar or a square matrix.' );
elseif size( Q, 1 ) ~= 1 && size( Q, 1 ) ~= length( x ),
    error( 'Sizes are incompatible.' );
end
if nargin < 3,
    v = 0;
elseif ndims( v ) > 1 || all( size( v ) ~= 1 ),
    error( 'The third argument must be a vector.' );
elseif size( Q, 1 ) ~= length( v ),
    error( 'Sizes are incompatible.' );
end
if nargin < 4,
    w = 0;
elseif numel( w ) > 1,
    error( 'The fourth argument must be a scalar.' );
end

x = x( : );
v = v( : );
y = real( x' * ( Q * x ) ) + v' * x + w;

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
