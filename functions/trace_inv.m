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

persistent P
if isempty( P ),
    P.nargs     = 1;
    P.args      = [];
    P.empty     = 0;
    P.constant  = @trace_inv_diag;
    P.diagonal  = @trace_inv_diag;
    P.affine    = @trace_inv_aff;
    P.structure = 'psdeig';
end
y = cvx_matrix_op( P, varargin );

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
