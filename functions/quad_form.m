function y = quad_form( x, Q )

%QUAD_FORM quadratic form.
%   QUAD_FORM(x,Q) is x'*((Q+Q'/2))*x = real(x'*Q*x).
%
%   x must be a row or column vector, and Q must either be a scalar or
%   a square matrix with the same number of rows as x(:).
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
%       If x is constant, then QUAD_FORM is affine in Q. The monotonicity
%       of QUAD_FORM depends on the precise values of x in this case, which
%       in turn govern whether the elements of Q can be convex, concave, or
%       affine. An error message will be generated if an inappropriate sum
%       of convex and concave terms occurs.

%
% Check sizes and types
%

error( nargchk( 2, 2, nargin ) );
sx = size( x );
if length( sx ) ~= 2 | all( sx > 1 ),
    error( 'The first argument must be a row or column.' );
elseif ndims( Q ) > 2 | size( Q, 1 ) ~= size( Q, 2 ),
    error( 'The second argument must be a scalar or a square matrix.' );
elseif size( Q, 1 ) ~= 1 & size( Q, 1 ) ~= length( x ),
    error( 'Sizes are incompatible.' );
end

x = x( : );
y = real( x' * ( Q * x ) );

% Copyright 2008 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
