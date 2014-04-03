function y = square( x, mode )

% SQUARE   Internal cvx version.
% Also implements SQUARE_POS and SQUARE_ABS.

persistent remap
if isempty( remap ),
    remap = cvx_remap( 'constant' );
    remap = remap + 2 * ( cvx_remap( 'l-valid' ) & ~remap );
    remap = remap + 3 * ( cvx_remap( 'r-affine' ) & ~remap );
    remap = remap + 4 * ( cvx_remap( 'p-convex' ) & ~remap );
    remap = remap + 5 * ( cvx_remap( 'n-concave' ) & ~remap );
end
if nargin == 2,
    switch mode,
        case 'abs', x = abs(x);
        case 'pos', x = pos(x);
    end
end
v = remap( cvx_classify( x ) );

%
% Perform the computations for each expression type separately
%

yt = [];
vu = sort( v(:) );
vu = vu([true;diff(vu)~=0]);
nv = length( vu );
nd = ndims( x ) + 1;
if nv ~= 1,
    y = cvx( size( x ), [] );
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
    sx = size( xt );

    %
    % Perform the computations
    %

    switch vk,
        case 0,
            % Invalid
            if nargin < 2, mode = ''; elseif isempty( mode ), mode = ''; else mode = [ '_', mode ]; end %#ok
            error( 'Disciplined convex programming error:\n    Illegal operation: square%s( {%s} ).', mode, cvx_class( xt, true, true ) );
        case 1,
            % Constant
            yt = cvx( cvx_constant( xt ) .^ 2 );
        case 2,
            % Monomial, posynomial
            yt = exp( 2 * log( xt ) );
        case {3,4,5},
            % Real affine, p-convex, n-concave
            cvx_begin
                epigraph variable yt( sx )
                if vk ~= 3, xt = cvx_accept_cvxccv( xt ); end
                { xt, 1, yt } == rotated_lorentz( sx, nd, 0 ); %#ok
                cvx_setnneg(yt);
            cvx_end
        otherwise,
            error( 'Shouldn''t be here.' );
    end

    %
    % Store the results
    %

    if nv == 1,
        y = yt;
    else
        y = cvx_subsasgn( y, t, yt );
    end

end
    
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
