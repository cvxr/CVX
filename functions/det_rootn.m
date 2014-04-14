function y = det_rootn( X )

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
    params.funcs  = { @det_rootn_cnst, @det_rootn_aff };
    params.square = true;
    params.name   = 'det_rootn';
end

try
    y = matrix_op( params, X );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function y = det_rootn_cnst( X )
[ psd, X, Z ] = cvx_check_psd( X, 'chol' ); %#ok
if ~psd,
	y = -Inf;
else
	y = geo_mean( diag( Z ) ) .^ 2;
end

function cvx_optval = det_rootn_aff( X )
n = size( X, 1 );
if nnz( X ) <= n && nnz( diag( X ) ) == nnz( X ),
    cvx_optval = geo_mean( diag( X ) );
else
    cvx_begin sdp
    	if isreal( X ),
	        variable Z(n,n) lower_triangular
	    else
	        variable Z(n,n) lower_triangular complex
	    end
        D = real( diag( Z ) );
        maximize( geo_mean( D ) );
        subject to
            [ diag( D ), Z' ; Z, X ] >= 0;
    cvx_end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
