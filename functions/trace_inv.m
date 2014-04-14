function y = trace_inv( X )

% TRACE_INV   Trace of the inverse of a PSD matrix.
%     For square matrix X, TRACE_INV(X) is TRACE(INV(X)) if X is Hermitian
%     or symmetric and positive definite; and +Inf otherwise. 
%
%     An error results if X is not a square matrix.
%
%     Disciplined convex programming information:
%         TRACE_INV is convex and nonmonotonic (at least with respect to
%         elementwise comparison), so its argument must be affine.

persistent params
if isempty( params ),
    params.funcs  = { @trace_inv_cnst, @trace_inv_real, @trace_inv_cplx };
    params.square = true;
    params.name   = 'trace_inv';
end

try
    y = matrix_op( params, X );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function z = trace_inv_cnst( X )
err = X - X';
X   = 0.5 * ( X + X' );
if norm( err, 'fro' ) > 8 * eps * norm( X, 'fro' ),
	z = Inf;
else
	z = eig( full( X ) );
	if any( z <= 0 ),
		z = Inf;
	else
		z = sum( 1.0 ./ z );
	end
end

function cvx_optval = trace_inv_real( X ) %#ok
sx = size(X);
cvx_begin sdp
    variable Y(sx) symmetric
    cvx_setnneg(diag(Y));
    minimize(trace(Y));
    [Y,eye(sx);eye(sx),X] >= 0; %#ok
cvx_end

function cvx_optval = trace_inv_cplx( X ) %#ok
sx = size(X);
cvx_begin sdp
    variable Y(sx) Hermitian
    cvx_setnneg(diag(Y));
    minimize(trace(Y));
    [Y,eye(sx);eye(sx),X] >= 0; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
