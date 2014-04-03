function y = pdom( x )

% PDOM(X) internal CVX implementation.

persistent remap
if isempty( remap ),
    remap = cvx_remap( 'nonnegative', 'p-nonconst' );
    remap = remap + 2 * cvx_remap( 'negative' );
    remap = remap + 3 * ( cvx_remap( 'n-concave' ) & ~remap );
    remap = remap + 4 * ( cvx_remap( 'r-affine' ) & ~remap );
    remap = remap + 5 * ( cvx_remap( 'concave' ) & ~remap );
end
v = remap( cvx_classify( x ) );

%
% Process each type of expression one piece at a time
%

yt = [];
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
    st = size( xt ); %#ok

    %
    % Perform the computations
    %

    switch vk,
        case 0,
            % Invalid
            error( 'Disciplined convex programming error:\n    Illegal operation: abs( {%s} ).', cvx_class( xt ) );
        case 1,
            % Nonnegative
            yt = xt;
        case 2,
            % Negative
            error( 'Disciplined convex programming error:\n    Trivially infeasible: p( %g ).', min( xt(:) ) );
        case 3,
            % Nonpositive
            yt = 0;
            warning( 'Disciplined convex programming error:\n    Almost certainly infeasible: p( {%s} ).', cvx_class( xt ) );
            cvx_begin
                xt >= 0; %#ok
            cvx_end
        case 4,
            % Real affine
            cvx_begin
                variable yt( st ) nonnegative
                yt == xt; %#ok
            cvx_end
        case 5,
            % Concave
            cvx_begin
                hypograph variable yt( st ) nonnegative
                yt <= xt; %#ok
            cvx_end
        case 6,
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