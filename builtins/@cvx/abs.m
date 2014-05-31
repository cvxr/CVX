function y = abs( x )

%Disciplined convex/geometric programming information for ABS:
%   ABS(X) is convex and nonmonotonic in X. Therefore, according to
%   the strict rules of DCP, X must be affine. However, because of
%   its special structure, CVX considers the sign of X as well. So,
%   for instance, if X is known to be nonnegative, then ABS(X)=X.

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'constant' }, { 'p_nonconst' }, { 'n_nonconst' }, ...
        { 'r_affine' }, { 'c_affine' } );
    P.funcs = { @abs_cnst, @abs_posn, @abs_negn, @abs_affn, @abs_affn };
end
y = cvx_unary_op( P, x );

function y = abs_cnst( x )
% Constant
y = builtin( 'abs', x );

function y = abs_posn( x )
% Positive any
y = x;

function y = abs_negn( x )
% Negative any
y = -x;

function y = abs_affn( x ) %#ok
% Affine
cvx_begin
    epigraph variable y( size(x) ) nonnegative_
    { x, y } == lorentz( size(x), 0, ~isreal(x) ); %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
