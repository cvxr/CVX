function z = matrix_frac( varargin )

%MATRIX_FRAC   Matrix fractional function.
%     MATRIX_FRAC(X,Y), where Y is a square matrix and X is a vector of the
%     same size, computes X'*inv(Y)*X if Y is Hermitian positive definite
%     and +Inf otherwise.
%
%     If X is a matrix with the same number of rows as Y, then MATRIX_FRAC
%     computes TRACE(X'*INV(Y)*X) if Y is Hermtian positive definite and
%     +Inf otherwise.
%
%     An error results if Y is not a square matrix, or the size of
%     the vector x does not match the size of matrix Y.
%
%     Disciplined convex programming information:
%         MATRIX_FRAC is convex and nonmonotonic (at least with respect to
%         elementwise comparison), so its argument must be affine.

persistent params
if isempty( params ),
    params.nargs     = 2;
    params.args      = [];
    params.empty     = 0;
    params.args      = @matrix_frac_args;
    params.constant  = @matrix_frac_cnst;
    params.diagonal  = [];
    params.affine    = @matrix_frac_aff;
    params.structure = 'chol';
end

try
    z = cvx_matrix_op( params, varargin );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function [ Y, X ] = matrix_frac_args( X, Y )
if size( X, 1 ) ~= size( Y, 1 ),
    cvx_throw( 'The number of rows in X and Y must match.' );
end

function z = matrix_frac_cnst( Z, X )
z = norm( Z' \ X, 'fro' ) .^ 2;

function cvx_optval = matrix_frac_aff( Y, X )
[m,n] = size(X); %#ok
if m == 1,
    cvx_optval = quad_over_lin( X, Y );
else
    cvx_begin sdp
        Q = [ Y, X ];
        variable Z(m+n,m+n) hermitian_if(Q) semidefinite
        minimize(trace(Z(m+1:end,m+1:end))) %#ok
        Z(1:m,:) == Q; %#ok
    cvx_end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
