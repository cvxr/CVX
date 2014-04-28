function y = det_rootn( varargin )

%DET_ROOTN nth-root of the determinant of an SPD matrix.
%   For a square matrix X, DET_ROOTN(X) returns
%       POW(DET(X),1/(size(X,1))
%   if X is symmetric (real) or Hermitian (complex) and positive semidefinite,
%   and -Inf otherwise.
%
%   This function can be used in many convex optimization problems that call for
%   LOG(DET(X)) instead. For example, if the objective function contains nothing
%   but LOG(DET(X)), it can be replaced with DET_ROOTN(X), and the same optimal 
%   point will be produced.
%
%   Disciplined convex programming information:
%       DET_ROOTN is concave and nonmonotonic; therefore, when used in
%       CVX specifications, its argument must be affine.

persistent params
if isempty( params ),
    params.nargs     = 1;
    params.args      = [];
    params.empty     = 1;
    params.constant  = @det_rootn_diag;
    params.diagonal  = @det_rootn_diag;
    params.affine    = @det_rootn_aff;
    params.structure = 'psdeig';
end

try
    y = cvx_matrix_op( params, varargin );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function y = det_rootn_diag( D )
y = geo_mean( D ) .^ 2;

function cvx_optval = det_rootn_aff( X ) %#ok
cvx_begin sdp
    hypograph variable z nonnegative_
    variable Z(size(X)) lower_triangular complex_if(X)
    D = real( diag( Z ) );
    geo_mean( D ) >= z;
    [ diag( D ), Z' ; Z, X ] >= 0;
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
