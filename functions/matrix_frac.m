function z = matrix_frac( X, Y )

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
    params.funcs  = { @matrix_frac_cnst, @matrix_frac_aff };
    params.square = [false,true];
    params.name   = 'trace_sqrtm';
end

try
	if size( X, 1 ) ~= size( Y, 1 ),
    	error( 'The number of rows in X and Y must match.' );
    end
    z = matrix_op( params, X, Y );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function z = matrix_frac_cnst( X, Y )
[ psd, Y, Z ] = cvx_check_psd( Y, 'chol' ); %#ok
if psd, z = norm( Z' \ X, 'fro' ) .^ 2;
else z = Inf; end

function cvx_optval = matrix_frac_aff( X, Y ) %#ok
n = size(X,1); %#ok
cvx_begin sdp
    variable z(n) assert_nonnegative
    minimize(sum(z))
    [ Y, X ; X', diag(z) ] >= 0; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
