function y = det_inv( varargin )

%DET_INV determinant of the inverse of an SPD matrix.
%   For a square matrix X, DET_INV(X) returns 1.0./DET(X) if X is symmetric
%   (real) or Hermitian (complex) and positive defininte, and +Inf otherwise.
%
%   This function can be used in many convex optimization problems that call
%   for LOG(DET(X)) instead. For example, if the objective function is
%      maximize(logdet(X))
%   then it can be replaced with
%      maximize(-det_inv(X))
%   and the same optimal point will be produced.
%
%   DET_INV(X,p) computes DET_INV(X)^p. p must be a positive real scalar.
%
%   Disciplined convex programming information:
%       DET_INV(X) is convex and nonmonotonic in X; therefore, when used in
%       CVX specifications, its argument must be affine.

persistent params
if isempty( params ),
    params.nargs     = 2;
    params.args      = @det_inv_args;
    params.empty     = 1;
    params.constant  = @det_inv_diag;
    params.diag      = @det_inv_diag;
    params.affine    = @det_inv_aff;
    params.structure = 'psdeig';
end

try
    y = cvx_matrix_op( params, varargin );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function [ X, p ] = det_inv_args( X, p )
if isempty( p ),
    p = 1;
elseif ~( isnumeric(p) && isreal(p) && numel(p)==1 && p>=0 ),
    error( 'CVX:ArgError', 'Second argument must be a positive scalar.' );
end

function y = det_inv_diag( D, p )
y = prod_inv( D, p );

function y = det_inv_aff( X, p )
cvx_begin sdp
    epigraph variable z nonnegative_
    variable Y(n,n) lower_triangular complex_if(X)
    prod_inv( real(diag(Y)), 2*p ) <= z;
    [ diag( D ), Y' ; Y, X ] >= 0;
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
