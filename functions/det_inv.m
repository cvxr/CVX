function y = det_inv( X, p )

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
    params.funcs  = { @det_inv_cnst, @det_inv_aff };
    params.square = true;
    params.name   = 'det_inv';
end

try
	if nargin < 2,
	    p = 1;
	elseif ~isnumeric( p ) || ~isreal( p ) || numel( p ) ~=  1 || p <= 0,
	    error( 'Second argument must be a positive scalar.' );
	end
    y = matrix_op( params, X, p );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function y = det_inv_cnst( X, p )
[psd,X,Z] = check_psd( X, 'chol' );	
if ~psd,
	y = Inf;
else
	y = prod( diag(Z) .^ ( -2 * p ) );
end

function y = det_inv_aff( X, p )
if nnz( X ) <= n && nnz( diag( X ) ) == nnz( X ),
	y = prod_inv( diag(X), p );
else
    cvx_begin sdp
        if isreal( X ),
	        variable Z(n,n) lower_triangular
        else
	        variable Z(n,n) lower_triangular complex
	    end
	    minimize(prod_inv(real(diag(Z)),2*p))
        [ diag( D ), Z' ; Z, X ] >= 0;
    cvx_end
end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
