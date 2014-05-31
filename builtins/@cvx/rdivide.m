function z = rdivide( x, y, oper )

%   Disciplined convex programming information for RDIVIDE:
%      For DCP purposes, the RDIVIDE division operator X./Y is identical
%      to X.*(1./Y). The right-hand term must be constant and non-zero;
%      and if the left-hand term is nonlinear, the constat must be real.
%   
%   Disciplined geometric programming information for RDIVIDE:
%      Terms in a left divide must have opposite log-curvature, so the
%      following products are permitted:
%         {log-convex} ./ {log-concave}  {log-concave} ./ {log-convex}
%         {log-affine} ./ {log-affine}
%      Note that log-affine expressions are both log-convex and 
%      log-concave.
%   
%   For vectors, matrices, and arrays, these rules are verified 
%   indepdently for each element.

persistent P
if isempty( P ),
    P.map = cvx_remap( ...
        ... % Invalid combinations: division by zero, division by posynomial,
        ... % negative of log-concave / gp monomial/posynomial
        { { 'valid' }, { 'zero', 'g_lconvex' } }, ...
        { { 'l_concave_', 'g_valid' }, { 'negative' } }, ...
        { { 'negative' }, { 'l_convex_', 'g_valid' } }, ...
        ... % Constant
        { { 'affine' }, { 'constant' } }, ...
        { { 'valid' }, { 'real' } }, ...
        ... % Geometric
        { { 'l_convex' },   { 'l_concave' } }, ...
        { { 'l_concave' },  { 'l_convex' }  }, ...
        ... % Implied inv_pos
        { { 'real' }, { 'p_concave' } }, ...
        { { 'real' }, { 'n_convex' }  }, ...
        [ 0, 0, 0, 1, 1, 2, 2, 3, 3 ] );
    P.funcs = { @rdivide_c, @rdivide_g, @rdivide_n };
    P.constant = [];
end
if nargin < 3, oper = './'; end
P.name = oper;
z = cvx_binary_op( P, x, y );

function z = rdivide_c( x, y )
% something ./ constant
z = bsxfun( @rdivide, cvx_basis( x ), cvx_constant( y ).' );
z = cvx( size(z,2), z );

function z = rdivide_g( x, y )
% geom ./ geom
z = exp( minus_nc( log( x ), log( y ) ) );

function z = rdivide_n( x, y )
% constant ./ something
% assumption: size(cvx_basis(x),1) == 1
z = bsxfun( @times, cvx_basis( x ), cvx_basis( recip( y ) ) );
z = cvx( size(z,2), z );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
