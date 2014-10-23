function y = norm_nuc( varargin )

%NORM_NUC   Nuclear norm of a matrix.
%   NORM_NUC(X) = SUM(SVD(X)). X must be a 2-D matrix, real or complex.
%
%   Disciplined convex programming information:
%       NORM_NUC(X) is convex and nonmontonic in X, so X must be affine.

persistent params
if isempty( params ),
	params.nargs     = 1;
	params.args      = [];
	params.empty     = 0;
	params.constant  = @norm_nuc_diag;
	params.diagonal  = @norm_nuc_diag;
	params.affine    = @norm_nuc_aff;
    params.structure = 'svd';
    params.name      = 'norm_nuc';
end

try
    y = cvx_matrix_op( params, varargin );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end        

function y = norm_nuc_diag( D )
y = sum(abs(D));

function cvx_optval = norm_nuc_aff( X ) %#ok
[ m, n ] = size( X ); %#ok
cvx_begin
	variable W(m+n,m+n) hermitian_if(X) semidefinite
	minimize(0.5*trace(W))
	W(1:m,m+1:end) == X;
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
