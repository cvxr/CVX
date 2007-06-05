function y = cvx_recip( x )
error( nargchk( 1, 1, nargin ) );

%RECIP    Reciprocal.
%
%   RECIP(X) is the elementwise inverse of X; i.e., 1.0 ./ X.
%       Unlike the native MATLAB division, division by zero is NOT
%       permitted and will result in an error.
%
%   Disciplined quadratic programming information:
%       When used in CVX expressions, X must be constant or log-affine.

%
% Determine the expression types
%

persistent remap
if isempty( remap ),
    remap_1 = cvx_remap( 'zero' );
    remap_2 = cvx_remap( 'constant' );
    remap_3 = cvx_remap( 'log-valid' );
    remap = remap_1 + 2 * ( remap_2 & ~remap_1 ) + 3 * ( remap_3 & ~remap_2 );
end
vr = remap( cvx_classify( x ) );
vu = unique( vr );
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
            error( sprintf( 'Disciplined convex programming error:\n    Cannot perform the operation recip( {%s} )', cvx_class( xt ) ) );
        case 1,
            % Zero
            error( sprintf( 'Disciplined convex programming error:\n    Division by zero.' ) );
        case 2,
            % Non-zero constant
            yt = 1.0 ./ xt;
        case 3,
            % Monomial, posynomial
            yt = exp( -log( xt ) );
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

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
