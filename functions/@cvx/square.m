function y = square( x )

%SQUARE    Square.
%   Internal CVX version. (Split from the constant version in case there is
%   a conflict with, say, the signal processing toolbox.)

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'real' }, { 'l_valid' }, ...
        { 'r_affine', 'p_convex', 'n_concave' } );
    P.funcs = { @square_cnst, @square_logv, @square_affn };
end
y = cvx_unary_op( P, x );

function y = square_cnst( x )
y = builtin( 'power', x, 2 );

function y = square_logv( x )
y = exp( 2 * log( x ) );

function y = square_affn( x ) %#ok
cvx_begin
    epigraph variable y( size(x) ) nonnegative_
    { cvx_linearize( x, 'abs' ), 0.5, y } == rotated_lorentz( size(x), 0 ); %#ok
cvx_end
    
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
