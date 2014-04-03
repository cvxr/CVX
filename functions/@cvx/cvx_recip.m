function y = cvx_recip( x )

%RECIP   Internal cvx version.

%
% Determine the expression types
%

error( nargchk( 1, 1, nargin ) ); %#ok
persistent remap
if isempty( remap ),
    remap = cvx_remap( 'constant' ) & ~cvx_remap( 'zero' );
    remap = remap + 2 * ( cvx_remap( 'l-valid' ) & ~remap );
    remap = remap + 3 * ( cvx_remap( 'p-concave' ) & ~remap );
    remap = remap + 4 * ( cvx_remap( 'n-convex' ) & ~remap );
end
vr = remap( cvx_classify( x ) );
vu = sort( vr(:) );
vu = vu([true;diff(vu)~=0]);
nv = length( vu );

%
% Process each result type one at a time
%

if nv ~= 1,
    y = cvx( x.size_, [] );
end
for k = 1 : nv,

    %
    % Select the category of expression to compute
    %

    if nv == 1,
        xt = x;
    else
        t = vr == vu( k );
        xt = cvx_subsref( x, t );
    end

    %
    % Perform the computations
    %

    switch vu( k ),
        case 0,
            % Invalid
            error( 'Disciplined convex programming error:\n    Cannot perform the operation recip( {%s} )', cvx_class( x, false, false, true ) );
        case 1,
            % Non-zero constant
            yt = cvx( 1.0 ./ cvx_constant( xt ) );
        case 2,
            % Monomial, posynomial
            yt = exp( -log( xt ) );
        case 3,
            % Positive concave
            yt = power( xt, -1 );
        case 4,
            % Negative convex
            yt = - power( - xt, -1 );
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
