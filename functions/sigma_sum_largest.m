function y = sigma_sum_largest( varargin )

% SIGMA_SUM_LARGEST   Sum of the k largest eigenvalues of a matrix.
%
%     For matrix X, SIGMA_SUM_LARGEST(X,K) is SUM_LARGEST(SVD(X),k).
%
%     Disciplined convex programming information:
%         SIGMA_SUM_LARGEST is convex and nonmonotonic (at least with 
%         respect to elementwise comparison), so its argument must be affine.

persistent params
if isempty( params ),
    params.nargs     = 2;
    params.args      = @sigma_sum_largest_args;
    params.empty     = 0;
    params.constant  = @sigma_sum_largest_cnst;
    params.diagonal  = @sigma_sum_largest_diag;
    params.affine    = @sigma_sum_largest_aff;
    params.structure = 'svd';
    params.name      = 'sigma_sum_largest';
end

try
    y = cvx_matrix_op( params, varargin );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function [ X, k ] = sigma_sum_largest_args( X, k )
if ~( isnumeric(k) && numel(k) == 1 && isreal(k) ),
    cvx_throw( 'Second input must be a constant real scalar.' );
end

function y = sigma_sum_largest_diag( D, k )
y = sum_largest( abs(D), k );

function cvx_optval = sigma_sum_largest_aff( X, k )
if k <= 0,
    y = 0;
elseif k <= 1,
    y = k * sigma_max( X );
elseif k >= min(size(X,1),size(X,2)),
    y = norm_nuc( X );
else
    [m,n] = size(X);
    cvx_begin sdp
        variable z nonnegative_
        variable S(m+n,m+n) hermitian_if(X) semidefinite
        z * eye(m+n) + S >= [zeros(m,m),X;X',zeros(n,n)]; %#ok
        minimize( k * z + trace( S ) );
    cvx_end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
