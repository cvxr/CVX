function y = trace_inv( varargin )

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
    params.nargs     = 1;
    params.args      = [];
    params.empty     = 0;
    params.constant  = @trace_inv_diag;
    params.diagonal  = @trace_inv_diag;
    params.affine    = @trace_inv_aff;
    params.structure = 'psdeig';
end

try
    y = cvx_matrix_op( params, varargin );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function z = trace_inv_diag( D )
z = sum( 1.0 ./ D );

function z = trace_inv_aff( X ) %#ok
cvx_begin sdp
    epigraph variable z nonnegative_
    variable Y(size(X)) hermitian_if(X)
    real(trace(Y)) <= z;
    [Y,eye(size(X));eye(size(X)),X] >= 0; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
