function y = trace_sqrtm( X )

% TRACE_SQRTM   Trace of the square root of a PSD matrix.
%     For square matrix X, TRACE_SQRTM(X) is TRACE(SQRTM(X)) if X is Hermitian
%     or symmetric and positive definite; and +Inf otherwise.
%
%     An error results if X is not a square matrix.
%
%     Disciplined convex programming information:
%         TRACE_SQRTM is convex and nonmonotonic (at least with respect to
%         elementwise comparison), so its argument must be affine.

persistent params
if isempty( params ),
    params.funcs  = { @trace_sqrtm_cnst, @trace_sqrtm_real, @trace_sqrtm_cplx };
    params.square = true;
    params.name   = 'trace_sqrtm';
end

try
    y = matrix_op( params, X );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function z = trace_sqrtm_cnst( X )
err = X - X';
X   = 0.5 * ( X + X' );
if norm( err, 'fro' ) > 8 * eps * norm( X, 'fro' ),
	z = Inf;
else
	z = eig( full( X ) );
	if any( z < 0 ),
		z = Inf;
	else
		z = sum( sqrt( z ) );
	end
end

function cvx_optval = trace_sqrtm_real( X ) %#ok
sx = size(X);
cvx_begin sdp
    variable Y(sx)
    cvx_setnneg(diag(Y));
    maximize(trace(Y))
    [ eye(sx), Y ; Y', X ] >= 0; %#ok
cvx_end

function cvx_optval = trace_sqrtm_cplx( X ) %#ok
sx = size(X);
cvx_begin sdp
    variable Y(sx) complex
    cvx_setnneg(real(diag(Y)));
    maximize(real(trace(Y)));
    [ eye(sx), Y ; Y', X ] >= 0; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
