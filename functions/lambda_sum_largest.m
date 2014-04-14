function y = lambda_sum_largest( X, k )

% LAMBDA_SUM_LARGEST   Sum of the k largest eigenvalues of a symmetric matrix.
%
%     For square matrix X, LAMBDA_SUM_LARGEST(X,K) is SUM_LARGEST(EIG(X),k)
%     if X is Hermitian or symmetric and real; and +Inf otherwise.
%
%     An error results if X is not a square matrix.
%
%     Disciplined convex programming information:
%         LAMBDA_SUM_LARGEST is convex and nonmonotonic (at least with 
%         respect to elementwise comparison), so its argument must be affine.

persistent params
if isempty( params ),
    params.funcs  = { @lambda_sum_largest_cnst, @lambda_sum_largest_aff };
    params.square = true;
    params.name   = 'lambda_sum_largest';
end

try
	if ~isnumeric(k) || numel(k) ~= 1 || ~isreal(k),
	    error( 'CVX:ArgError', 'Second input must be a constant real scalar.' );
	end
    y = matrix_op( params, X );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function y = lambda_sum_largest_cnst( X, k )
[ psd, X ] = cvx_check_psd( X, 'sym' );
if ~psd,
	y = Inf;
elseif k <= 0,
	y = 0;
elseif k >= size(X,1),
	y = trace(X);
else
	y = eig( full( X ) );
	y = sum(y(1:floor(k))) + (k-floor(k)) * y(ceil(k));
end

function cvx_optval = lambda_sum_largest_aff( X, k )
n = size(X,1);
if k <= 0,
	cvx_begin
		X == X'; %#ok
	cvx_end	
elseif k >= n,
	cvx_begin
		X == X'; %#ok
	cvx_end
	cvx_optval = trace( X );
else
    cvx_begin sdp
    	variable z
    	if isreal(X),
    		variable S(n,n) symmetric
    	else
    		variable S(n,n) hermitian
    	end
    	S >= 0; %#ok
    	z * eye(n) + S >= X; %#ok
        minimize( k * z + trace( S ) );
    cvx_end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.