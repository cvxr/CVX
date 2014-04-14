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

persistent mparam rparam lparam
if isempty( mparam ),
    mparam.map = max( 0, cvx_remap( ...
        { { 'negative' },   { 'l_concave_' } }, ... % negative of log-concave
        { { 'l_concave_' }, { 'negative' }   }, ... % negative of log-concave
        { { 'constant' },   { 'affine' }     }, ... % constant * affine/constant
        { { 'real' },       { 'valid'  }     }, ... % real * valid
        { { 'affine' },     { 'constant' }   }, ... % affine * constant
        { { 'valid' },      { 'real'  }      }, ... % valid * real
        { { 'affine' }                       }, ... % potential quadratic form
        { { 'p_convex', 'n_concave' }        }, ... % potential quadratic form
        { { 'l_convex' }                     }, ... % log-convex
        { { 'l_concave' }                    }, ... % log-concave
        [ -1, -1, 1, 1, 2, 2, 3, 3, 4, 4 ] ) );
    rparam.map = max( 0, cvx_remap( ...
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
        [ -1, -1, -1, 1, 1, 1, 2, 2, 4, 4 ] ) );
    lparam.map = rparam.map';
    mparam.funcs = { @times_1m, @times_2m, @times_3, @times_4m };
    rparam.funcs = { @times_1l, @times_2l, [], @times_4l };
    lparam.funcs = { @times_1r, @times_2r, [], @times_4r };
end

try
    if nargin < 3,
        oper = '.*';
    end
    switch oper,
        case { '*', '.*' },
            mparam.name = oper;
            param = mparam;
        case { '/', './' },
            rparam.name = oper;
            param = rparam;
        case { '\', '.\' },
            lparam.name = oper;
            param = lparam;
    end
    z = binary_op( param, x, y );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end

function z = times_1m( x, y )
% constant .* something
if isa( x, 'cvx' ), x = x.basis_;
else x = x'; end
z = bsxfun( @times, x, y.basis_ );
z = cvx( [size(z,2),1], z );

function z = times_1r( x, y )
% constant ./ something
if isa( x, 'cvx' ), x = x.basis_;
else x = x'; end
y = recip( y );
z = bsxfun( @times, x', y.basis_ );
z = cvx( [size(z,2),1], z );

function z = times_1l( x, y )
% something .\ constant
x = recip( x );
if isa( y, 'cvx' ), y = y.basis_;
else y = y'; end
z = bsxfun( @ldivide, x.basis_, y );
z = cvx( [size(z,2),1], z );

function z = times_2m( x, y )
% something .* constant
if isa( y, 'cvx' ), y = y.basis_;
else y = y'; end
z = bsxfun( @times, x.basis_, y );
z = cvx( [size(z,2),1], z );

function z = times_2r( x, y )
% something ./ constant
if isa( y, 'cvx' ), y = y.basis_;
else y = y'; end
z = bsxfun( @rdivide, x.basis_, y );
z = cvx( [size(z,2),1], z );

function z = times_2l( x, y )
% constant .\ something
if isa( x, 'cvx' ), x = x.basis_;
else x = x'; end
z = bsxfun( @ldivide, x, y.basis_ );
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

function z = times_4m( x, y )
% geom .* geom
z = exp( log( x ) + log( y ) );

function z = times_4r( x, y )
% geom ./ geom
z = exp( log( x ) - log( y ) );

function z = times_4l( x, y )
% geom .\ geom
z = exp( log( y ) - log( x ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
