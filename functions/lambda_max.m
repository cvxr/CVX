function y = lambda_max( X )

% LAMBDA_MAX    Maximum eigenvalue of a symmetric matrix.
%     For square matrix X, LAMBDA_MAX(X) is MAX(EIG(X)) if X is Hermitian
%     or symmetric and real; and +Inf otherwise. 
%
%     An error results if X is not a square matrix.
%
%     Disciplined convex programming information:
%         LAMBDA_MAX is convex and nonmonotonic (at least with respect to
%         elementwise comparison), so its argument must be affine.

persistent params
if isempty( params ),
    params.funcs  = { @lambda_max_cnst, @lambda_max_aff };
    params.square = true;
    params.name   = 'lambda_max';
end

try
    y = matrix_op( params, X );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw(exc);
    else rethrow(exc); end
end

function y = lambda_max_cnst( X )
[psd,X] = cvx_check_psd( X, 'sym' );
if psd,
    y = max(eig(full(X)));
else
    y = Inf; 
end

function z = lambda_max_aff( X ) %#ok
n = size( X, 1 );
cvx_begin sdp
    epigraph variable z
    z * eye( n ) >= X; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
