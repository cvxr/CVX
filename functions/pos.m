function y = pos( x )

% POS    Positive part.
%    POS(X) = MAX(X,0). X must be real.
%
%    Disciplined convex programming information:
%        POS(X) is convex and nondecreasing in X. Thus when used in CVX
%        expressions, X must be convex (or affine).

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'real' }, { 'p_nonconst' }, { 'n_nonconst' }, ...
        { 'r_affine', 'convex' } );
    P.funcs = { @pos_real, @pos_posn, @pos_negn, @pos_cnvx };
end
y = cvx_unary_op( P, x );

function y = pos_real( x )
y = max( x, 0 );

function y = pos_posn( x )
y = x;

function y = pos_negn( x )
y = cvx( size(x), [] );

function y = pos_cnvx( x ) %#ok
cvx_begin
    epigraph variable y(size(x)) nonnegative
    x <= y; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
