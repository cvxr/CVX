function z = power( x, y )

%
% Check sizes
%

sx = size( x );
sy = size( y );
xs = all( sx == 1 );
ys = all( sy == 1 );
if xs,
    sz = sy;
elseif ys | isequal( sx, sy ),
    sz = sx;
else
    error( 'Matrix dimensions must agree.' );
end

%
% Determine the expression types
%

% Classifications:
% 1  - negative constant
% 2  - zero
% 3  - positive constant
% 4  - complex constant
% 5  - concave
% 6  - real affine
% 7  - convex
% 8  - complex affine
% 9  - log concave
% 10 - log affine
% 11 - log convex monomial
% 12 - log convex posynomial
% 13 - invalid

persistent remap
if isempty( remap ),
    % Constant .^ constant
    remap_1 = cvx_remap( 'constant' );
    remap_1 = remap_1' * remap_1;

    % Zero .^ convex
    remap_2 = cvx_remap( 'zero' )' * cvx_remap( 'convex' );

    % Positive .^ convex
    remap_3 = cvx_remap( 'positive' )' * cvx_remap( 'convex' );
    
    % Concave/affine ^ constant: potential quadratic form or root
    remap_4 = cvx_remap( 'concave' )' * cvx_remap( 'positive' );

    % log-concave/convex ^ zero
    remap_5 = cvx_remap( 'log-convex', 'log-concave' )' * cvx_remap( 'zero' );

    % log-concave/convex ^ nonzero
    remap_6 = cvx_remap( 'log-convex', 'log-concave' )' * cvx_remap( 'nonzero' );

    remap = remap_1 + ( 2 * remap_2 + 3 * remap_3 + 4 * remap_4 + 5 * remap_5 + 6 * remap_6 ) .* ~remap_1;
end
vx = cvx_classify( x );
vy = cvx_classify( y );
vr = remap( vx + size( remap, 1 ) * ( vy - 1 ) );
vu = unique( vr );
nv = length( vu );

%
% Perform the individual computations and combine
%

x  = cvx( x );
y  = cvx( y );
xt = x;
yt = y;
if nv ~= 1,
    z = cvx( sz, [] );
end
for k = 1 : nv,

    %
    % Select the category of expression to compute
    %

    if nv ~= 1,
        t = vr == vu( k );
        if ~xs, xt = cvx_subsref( x, t ); sz = size( xt ); end
        if ~ys, yt = cvx_subsref( y, t ); sz = size( yt ); end
    end

    %
    % The computational kernels
    %

    switch vu( k ),
        case 0,
            % Invalid
            error( sprintf( 'Disciplined convex programming error:\n    Cannot perform the operation {%s}.^{%s}', cvx_class( xt ), cvx_class( yt ) ) );
        case 1,
            % constant .^ constant
            cvx_optval = cvx_constant( xt ) .^ cvx_constant( yt );
        case 2,
            % zero .^ convex
            cvx_optval = zeros( sz );
        case 3,
            % positive .^ convex
            cvx_optval = exp( log( cvx_constant( xt ) ) .* yt );
        case 4,
            % concave/affine ^ positive constant --- potential quadratic form or root
            yt = cvx_constant( yt );
            t1 = ( yt > 0 ) & ( yt <= 1 );
            t2 = ( yt > 1 ) & ( rem( yt, 2 ) == 0 );
            t3 = ~t1 & ~t2;
            if nnz( t3 ),
                if ~cvx_isaffine( xt ),
                    error( sprintf( ...
                        [ 'Disciplined convex programming error:\n', ...
                          '    Cannot perform {concave}.^k unless 0<k<=1.' ] ) );
                else
                    error( sprintf( ...
                        [ 'Disciplined convex programming error:\n', ...
                          '    Cannot perform {affine}.^k unless 0<k<=1 or k is even.' ] ) );
                end
            elseif all( t1( : ) ),
                cvx_optval = pow_pos( xt, yt );
            elseif all( t2( : ) ),
                cvx_optval = pow_pos( square( xt ), 0.5 * yt );
            else
                if xs, xt = xt * ones( sz ); end
                xt = cvx_subsasgn( xt, square( cvx_subsref( xt, t2 ) ) );
                if ys, yt = yt * ones( sz ); end
                yt = cvx_subsasgn( yt, 0.5 * cvx_subsref( yt, t2 ) );
                cvx_optval = pow_pos( xt, yt );
            end
        case 5,
            % { log concave, log affine, log convex } ^ zero
            cvx_optval = ones( sz );
        case 6,
            % { log concave, log affine, log convex } ^ real
            cvx_optval = exp( log( xt ) .* cvx_constant( yt ) );
        otherwise,
            error( 'Shouldn''t be here.' );
    end

    %
    % Store the results
    %

    if nv == 1,
        z = cvx_optval;
    else
        z = cvx_subsasgn( z, t, cvx_optval );
    end

end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
