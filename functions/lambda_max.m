function y = lambda_max( varargin )

% LAMBDA_MAX    Maximum eigenvalue of a symmetric matrix.
%     For square matrix X, LAMBDA_MAX(X) is MAX(EIG(X)) if X is Hermitian
%     or symmetric and real; and +Inf otherwise. 
%
%     An error results if X is not a square matrix.
%
%     Disciplined convex programming information:
%         LAMBDA_MAX is convex and nonmonotonic (at least with respect to
%         elementwise comparison), so its argument must be affine.

persistent params
if isempty( params ),
	params.nargs     = 1;
	params.args      = [];
	params.empty     = [];
	params.constant  = @lambda_max_diag;
	params.diagonal  = @lambda_max_diag;
	params.structure = 'eig';
end

try
    y = cvx_matrix_op( params, varargin );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function y = lambda_max_diag( D )
y = max( D );	

function z = lambda_max_aff( X ) %#ok
cvx_begin sdp
    epigraph variable z
    z * eye(size(x)) >= X; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
