function y = trace_sqrtm( varargin )

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
    params.nargs     = 1;
    params.args      = [];
    params.empty     = 0;
    params.constant  = @trace_sqrtm_diag;
    params.diagonal  = @trace_sqrtm_diag;
    params.affine    = @trace_sqrtm_aff;
    params.structure = 'psdeig';
end

try
    y = cvx_matrix_op( params, varargin );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function z = trace_sqrtm_diag( D )
z = sum( sqrt(D) );

function z = trace_sqrtm_aff( X ) %#ok
cvx_begin sdp
    hypograph variable z nonnegative_
    variable Y(size(X)) complex_if(X)
    real(trace(Y)) >= z; %#ok
    [ eye(size(X)), Y ; Y', X ] >= 0; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
