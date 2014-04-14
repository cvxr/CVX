function y = norm_nuc( X )

%NORM_NUC   Nuclear norm of a matrix.
%   NORM_NUC(X) = SUM(SVD(X)). X must be a 2-D matrix, real or complex.
%
%   Disciplined convex programming information:
%       NORM_NUC(X) is convex and nonmontonic in X, so X must be affine.

persistent params
if isempty( params ),
    params.funcs  = { @norm_nu_cnst, @norm_nuc_real, @norm_nuc_cplx };
    params.square = false;
    params.name   = 'norm_nuc';
end

try
    y = matrix_op( params, X );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end        

function y = norm_nuc_cnst( X )
y = sum(svd(X));

function cvx_optval = norm_nuc_real( X ) %#ok
[ m, n ] = size( X ); %#ok
cvx_begin
	variable W(m+n,m+n) semidefinite
	minimize(0.5*trace(W))
	W(1:m,m+1:end) == X;
cvx_end

function cvx_optval = norm_nuc_cplx( X ) %#ok
[ m, n ] = size( X ); %#ok
cvx_begin
	variable W(m+n,m+n) hermitian semidefinite
	minimize(0.5*trace(W))
	W(1:m,m+1:m+n) == X;
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
