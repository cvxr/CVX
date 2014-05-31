function y = lambda_sum_largest( varargin )

% LAMBDA_SUM_LARGEST   Sum of the k largest eigenvalues of a symmetric matrix.
%
%     For square matrix X, LAMBDA_SUM_LARGEST(X,K) is SUM_LARGEST(EIG(X),k)
%     if X is Hermitian or symmetric and real; and +Inf otherwise.
%
%     An error results if X is not a square matrix.
%
%     Disciplined convex programming information:
%         LAMBDA_SUM_LARGEST is convex and nonmonotonic (at least with 
%         respect to elementwise comparison), so its argument must be affine.

persistent params
if isempty( params ),
    params.nargs     = 2;
    params.args      = @lambda_sum_largest_args;
    params.empty     = 0;
    params.constant  = @lambda_sum_largest_diag;
    params.diagonal  = @lambda_sum_largest_diag;
    params.structure = 'eig';
    params.name      = 'lambda_sum_largest';
end

try
    y = cvx_matrix_op( params, varargin );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function [ X, k ] = lambda_sum_largest_args( X, k )
if ~( isnumeric(k) && numel(k)==1 && isreal(k) ),
    cvx_throw( 'Second input must be a constant real scalar.' );
end

function y = lambda_sum_largest_diag( D, k )
y = sum_largest( D, k );

function y = lambda_sum_largest_aff( X, k )
if k <= 0,
    y = 0;
elseif k <= 1,
    y = k * lambda_max( X );
elseif k >= size(X,1),
    y = trace(X);
else
    cvx_begin sdp
    	variable z
        epigraph variable y
        variable S(size(X)) hermitian_if(X) semidefinite
    	z * eye(n) + S >= X; %#ok
        k * z + trace( S ) <= y; %#ok
    cvx_end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
