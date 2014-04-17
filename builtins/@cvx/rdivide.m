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

persistent params
if isempty( params ),
    params.map = cvx_remap( ...
        { { 'valid' },      { 'zero' }      }, ... % division by zero
        { { 'l_concave_' }, { 'negative' }  }, ... % negative of log-concave
        { { 'negative' },   { 'l_convex_' } }, ... % negative of 1/log-convex
        { { 'positive' },   { 'l_valid' }   }, ... % constant / non-constant
        { { 'real' },       { 'p_concave' } }, ... % constant / non-constant
        { { 'real' },       { 'n_convex' }  }, ... % constant / non-constant
        { { 'affine' },     { 'constant' }  }, ... % non-constant / constant
        { { 'valid' },      { 'real' }      }, ... % non-constant / constant
        { { 'l_convex' },   { 'l_concave' } }, ... % log-convex / log-concave
        { { 'l_concave' },  { 'l_convex' }  }, ... % log-concave / log-convex
        [ 0, 0, 0, 1, 1, 1, 2, 2, 3, 3 ] );
    params.funcs = { @rdivide_1, @rdivide_2, @rdivide_3 };
end

try
    if nargin < 3, oper = './'; end
    params.name = oper;
    z = cvx_binary_op( params, x, y );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end

function z = rdivide_1( x, y )
% constant ./ something
% assumption: size(cvx_basis(x),1) == 1
z = bsxfun( @times, cvx_basis( x ), cvx_basis( recip( y ) ) );
z = cvx( [size(z,2),1], z );

function z = rdivide_2( x, y )
% something ./ constant
% assumption: size(cvx_basis(y),1) == 1
z = bsxfun( @rdivide, cvx_basis( x ), cvx_basis( y ) );
z = cvx( [size(z,2),1], z );

function z = rdivide_3( x, y )
% geom ./ geom
z = exp( minus_nc( log( x ), log( y ) ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
