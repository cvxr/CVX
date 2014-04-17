function z = times( x, y, oper )

%   Disciplined convex programming information for TIMES:
%      In general, disciplined convex programs must obey the "no-product
%      rule" which states that two non-constant expressions cannot be 
%      multiplied together. Under this rule, the following combinations
%      are allowed:
%         {constant} .* {non-constant}  {non-constant} .* {constant}
%      Furthermore, if the non-constant expression is nonlinear, then 
%      the constant term must be real.
%   
%      A lone exception to the no-product rule is made for quadratic 
%      forms: two affine expressions may be multiplied together if the 
%      result is convex or concave. For example, the construction
%         variable x(n)
%         x.*x  <= 1;
%      would be permitted because each element of x.*x is convex.
%   
%   Disciplined geometric programming information for TIMES:
%      Both terms in a multiplication must have the same log-curvature, 
%      so the following products are permitted:
%         {log-convex} .* {log-convex}  {log-concave} .* {log-concave}
%         {log-affine} .* {log-affine}
%      Note that log-affine expressions are both log-convex and
%      log-concave.
%   
%   For vectors, matrices, and arrays, these rules are verified 
%   indepdently for each element.

persistent params
if isempty( params ),
    params.map = cvx_remap( ...
        { { 'negative' },   { 'l_concave_' } }, ... % negative of log-concave
        { { 'l_concave_' }, { 'negative' }   }, ... % negative of log-concave
        { { 'constant' },   { 'affine' }     }, ... % constant * affine/constant
        { { 'real' },       { 'valid'  }     }, ... % real * valid
        { { 'affine' },     { 'constant' }   }, ... % affine * constant
        { { 'valid' },      { 'real'  }      }, ... % valid * real
        { { 'affine' }                       }, ... % potential quadratic form
        { { 'p_convex' }                     }, ... % potential quadratic form
        { { 'l_convex' }                     }, ... % log-convex
        { { 'l_concave', 'concave' }         }, ... % log-concave
        [ 0, 0, 1, 1, 2, 2, 3, 3, 4, 4 ] );
    params.funcs = { @times_12, @times_12, @times_3, @times_4 };
end

try
    if nargin < 3, oper = '.*'; end
    params.name = oper;
    z = cvx_binary_op( params, x, y );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end

function z = times_12( x, y )
% constant .* something OR something .* constant
% assumption: size(cvx_basis(x),1) == 1 OR size(cvx_basis(y),1)==1
z = bsxfun( @times, cvx_basis( x ), cvx_basis( y ) );
z = cvx( [size(z,2),1], z );

function z = times_3( x, y )
% affine .* affine
% p_convex .* p_convex
% n_concave .* n_concave
xA = x.basis_; 
yA = y.basis_;
mm = max( size( xA, 1 ), size( yA, 1 ) );
xA( end + 1 : mm, : ) = 0;
yA( end + 1 : mm, : ) = 0;
xB  = xA( 2 : end, : );
yB  = yA( 2 : end, : );
cyB = conj( yB );
alpha = sum( yB .* cyB, 1 );
alpha = sum( bsxfun( @times, xB, yB ), 1 ) ./ max( alpha, realmin );
if nnz( xB - bsxfun( @times, alpha, cyB ) > 2 * eps * sqrt( conj( xB ) .* xB ) ),
    error( 'Disciplined convex programming error:\n    Invalid quadratic form(s): not a square.\n', 1 ); %#ok
elseif any( abs(imag(alpha)) > 2 * eps * abs(real(alpha)) ),
    error( 'Disciplined convex programming error:\n    Invalid quadratic form(s): not real.\n', 1 ); %#ok
else
    if isreal( xB ),
        z = square( y );
    else
        z = square( abs( y ) );
    end
    z.basis_ = bsxfun( @times, alpha, z.basis_ );
    offset = xA(1,:) - alpha .* conj( yA(1,:) );
    if any( offset ),
        z = z + offset.' .* y;
    end
end

function z = times_4( x, y )
% geom .* geom
z = exp( plus_nc( log( x ), log( y ) ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
