function y = sigma_max( X )

%SIGMA_MAX    Maximum singular value.
%   SIGMA_MAX(X) returns the maximum singular value of X. X must be a 2-D
%   matrix, real or complex. SIGMA_MAX(X) is synonymous with NORM(X).
%
%   Disciplined convex programming information:
%       SIGMA_MAX(X) is convex and nonmontonic in X, so X must be affine.

persistent params
if isempty( params ),
    params.funcs  = { @sigma_max_cnst, @sigma_max_aff };
    params.square = false;
    params.name   = 'sigma_max';
end

try
    y = matrix_op( params, X );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function y = sigma_max_cnst( X )
y = norm( X );

function z = sigma_max_aff( X )
[ m, n ] = size( X );
cvx_begin sdp
    epigraph variable z assert_nonnegative
    z * speye(m,n) >= [zeros(m,m),X;X',zeros(n,n)]; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
