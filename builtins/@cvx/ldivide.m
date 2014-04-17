function z = ldivide( x, y, oper )

%Disciplined convex programming information for LDIVIDE:
%   For DCP purposes, the LDIVIDE division operator X.\Y is equivalent to
%   (1./X).*Y. The left-hand term must be constant and non-zero; and if the
%   right-hand term is nonlinear, t constant must also be real.
%
%Disciplined geometric programming information for LDIVIDE:
%   Terms in a left divide must have opposite log-curvature, so the
%   following products are permitted:
%      {log-convex} .\ {log-concave}  {log-concave} .\ {log-convex}
%      {log-affine} .\ {log-affine}
%   Note that log-affine expressions are both log-convex and log-concave.
%
%For vectors, matrices, and arrays, these rules are verified indepdently
%for each element.

persistent params
if isempty( params ),
    params.map = cvx_remap( ...
        { { 'zero' },      { 'valid' }      }, ... % division by zero
        { { 'negative' },  { 'l_concave_' } }, ... % negative of log-concave
        { { 'l_convex_' }, { 'negative' }   }, ... % negative of recip( log-convex )
        { { 'l_valid' },   { 'positive' }   }, ... % log-valid    \ constant
        { { 'p_concave' }, { 'real' }       }, ... % p-concave    \ real
        { { 'n_convex' },  { 'real' }       }, ... % n-convex     \ real
        { { 'constant' },  { 'affine' }     }, ... % nonzero      \ affine
        { { 'real' },      { 'valid' }      }, ... % nonzero real \ anything
        { { 'l_concave' }, { 'l_convex' }   }, ... % log-concave  \ log-convex
        { { 'l_convex' },  { 'l_concave' }  }, ... % log-convex   \ log-concave
        [ 0, 0, 0, 1, 1, 1, 2, 2, 3, 3 ] );
    params.funcs = { @ldivide_1, @ldivide_2, @ldivide_3 };
end

try
    if nargin < 3, oper = '.\'; end
    params.name = oper;
    z = cvx_binary_op( params, x, y );
catch exc
	if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
	else rethrow( exc ); end
end

function z = ldivide_1( x, y )
% something .\ constant
% assumption: size(cvx_basis(y),1) == 1
z = bsxfun( @ldivide, cvx_basis( recip( x ) ), cvx_basis( y ) );
z = cvx( [size(z,2),1], z );

function z = ldivide_2( x, y )
% constant .\ something
% assumption: size(cvx_basis(x),1) == 1
z = bsxfun( @ldivide, cvx_basis( x ), cvx_basis( y ) );
z = cvx( [size(z,2),1], z );

function z = ldivide_3( x, y )
% geom .\ geom
z = exp( minus_nc( log( y ), log( x ) ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
