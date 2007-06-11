function y = inv_pos( x )
error( nargchk( 1, 1, nargin ) );

% INV_POS   Inverse of a positive quantity.
%     For a real scalar X, INV_POS(X) returns 1/X if X is positive, and
%     +Inf otherwise.
%
%     For matrices and N-D arrays, the function is applied to each element.
%
%     Disciplined convex programming information:
%         INV_POS is convex and nonincreasing; therefore, when used in CVX
%         specifications, its argument must be concave (or affine).

%
% Determine the expression typess
% 0 : convex, invalid
% 1 : complex
% 2 : real constant
% 3 : concave
% 4 : monomial / posynomial
%

persistent remap
if isempty( remap ),
    remap_1 = cvx_remap( 'complex' );
    remap_2 = cvx_remap( 'real' );
    remap_3 = cvx_remap( 'concave' ) & ~remap_2;
    remap_4 = cvx_remap( 'monomial', 'posynomial' );
    remap = remap_1 + 2 * remap_2 + 3 * remap_3 + 4 * remap_4;
end
v = remap( cvx_classify( x ) );

%
% Process each type of expression one piece at a times
%

vu = unique( v );
nv = length( vu );
sx = size( x );
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
    else,
        t = v == vk;
        xt = cvx_subsref( x, t );
    end
    
    %
    % Perform the computations
    %
    
    switch vk,
        case 0,
            % Invalid
            error( sprintf( 'Disciplined convex programming error:\n    Illegal operation: inv_pos( {%s} )', cvx_class( xt ) ) ); 
        case 1,
            % Complex
            error( 'The argument to inv_pos must be real.' );
        case 2,
            % Constant
            yt = inv_pos( cvx_constant( xt ) );
        case 3,
            % Concave
            cvx_begin
                epigraph variable yt( size( xt ) )
                quad_over_lin( 1, xt ) <= yt;
            cvx_end
            yt = cvx_optval;
        case 4,
            % Monomial / psosynomial
            yt = exp( -log( xt ) );
        otherwise,
            error( 'Shouldn''t be here.' );
    end
    
    %
    % Store the results
    %
    
    if nv == 1,
        y = yt;
    else,
        y = cvx_subsasgn( y, t, yt );
    end
    
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

