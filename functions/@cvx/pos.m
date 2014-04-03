function y = pos( x )

% POS(X) internal CVX implementation.

persistent remap
if isempty( remap ),
    remap = cvx_remap( 'real' );
    remap = remap + 2 * ( cvx_remap( 'p-nonconst' ) & ~remap );
    remap = remap + 3 * ( cvx_remap( 'n-nonconst' ) & ~remap );
    remap = remap + 4 * ( cvx_remap( 'r-affine', 'convex' ) & ~remap );
end
v = remap( cvx_classify( x ) );

%
% Process each type of expression one piece at a time
%

vu = sort( v(:) );
vu = vu([true;diff(vu)~=0]);
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
            error( 'Disciplined convex programming error:\n    Illegal operation: abs( {%s} ).', cvx_class( xt ) );
        case 1,
            % Constant
            yt = cvx( max( cvx_constant( xt ), 0 ) );
        case 2,
            % Positive any
            yt = xt;
        case 3,
            % Negative any
            yt = -xt;
        case 4,
            % Real affine, convex
            st = size( xt ); %#ok
            cvx_begin
                epigraph variable yt( st )
                xt <= yt; %#ok
                0  <= yt; %#ok
                cvx_setnneg( yt );
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
