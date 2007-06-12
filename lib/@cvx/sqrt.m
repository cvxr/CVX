function cvx_optval = sqrt( x )
error( nargchk( 1, 1, nargin ) );

%SQRT    Square root.
%
%   SQRT(X) is the square root of the elements of X. If X is numeric, then
%   complex results are returned for those elements of X which are complex
%   or real and negative. If X is a CVX expression, then X must be real,
%   the function effectively constrains X to be nonnegative.
%
%   Disciplined quadratic programming information:
%       SQRT(X) is concave and nondecreasing in X. Thus when used in CVX
%       expressions, X must be concave (or affine).

%
% Determine the expression types
%

% 0 : affine complex, convex, invalid
% 1 : constant
% 2 : concave, real affine
persistent remap
if isempty( remap ),
    remap_1 = cvx_remap( 'constant' );
    remap_2 = cvx_remap( 'real-affine' );
    remap_3 = cvx_remap( 'log-convex', 'log-concave' );
    remap = remap_1 + ( 2 * remap_2 + 3 * remap_3 ) .* ~remap_1;
end
v = remap( cvx_classify( x ) );

%
% Perform the computations for each expression type separately
%

vu = unique( v );
nv = length( vu );
sx = x.size_;
if nv ~= 1,
    y = cvx( sx, [] );
end
for k = 1 : nv,

    %
    % Select the category of expression to compute
    %

    vk = vu( k );
    if nv == 1,
        xt = x;
    else
        t = v == vk;
        xt = cvx_subsref( x, t );
    end

    %
    % Perform the computations
    %

    switch vk,
        case 0,
            % Invalid
            error( sprintf( 'Disciplined convex programming error:\n    Illegal operation: sqrt( {%s} ).', cvx_class( xt, true ) ) );
        case 1,
            % Constant
            cvx_optval = builtin( 'sqrt', cvx_constant( xt ) );
        case 2,
            % Real affine, concave
            st = size( xt );
            cvx_begin
                hypograph variable w( st );
                square( w ) <= xt;
            cvx_end
        case 3,
            % Monomial, posynomial
            cvx_optval = exp( 0.5 * log( xt ) );
        otherwise,
            error( 'Shouldn''t be here.' );
    end

    %
    % Store the results
    %

    if nv == 1,
        y = cvx_optval;
    else
        y = cvx_subsasgn( y, t, cvx_optval );
    end

end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
