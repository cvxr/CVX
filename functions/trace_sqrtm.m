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

persistent P
if isempty( P ),
    P.nargs     = 1;
    P.args      = [];
    P.empty     = 0;
    P.constant  = @trace_sqrtm_diag;
    P.diagonal  = @trace_sqrtm_diag;
    P.affine    = @trace_sqrtm_aff;
    P.structure = 'psdeig';
end
y = cvx_matrix_op( P, varargin );

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
