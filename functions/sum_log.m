function y = sum_log( varargin )

%SUM_LOG   Sum of logarithms.
%   For vectors, SUM_LOG(X) is the sum of the logarithms of the elements of
%   the vector; i.e., SUM(LOG(X)). If any of the elements of the vector are
%   nonnegative, then the result is -Inf.
%
%   For matrices, SUM_LOG(X) is a row vector containing the application of
%   SUM_LOG to each column. For N-D arrays, the SUM_LOG is applied to the
%   first non-singleton dimension of X.
%
%   SUM_LOG(X,DIM) takes the sum along the dimension DIM of X.
%
%   SUM_LOG(X) could also be written SUM(LOG(X)). However, this version
%   should be more efficient, because it involves only one logarithm.
%
%   Disciplined convex programming information:
%       SUM_LOG(X) is concave and nondecreasing in X. Therefore, when used
%       in CVX expressions, X must be concave. X must be real.

cvx_expert_check( 'sum_log', varargin{1} );

persistent params
if isempty( params ),
    params.map = cvx_remap( { 'real' ; 'l_convex' ; 'concave' ; 'l_concave' } );
    params.map = bsxfun( @and, params.map, ~cvx_remap( 'nonpositive' ) );
    params.funcs = { @sum_log_c, @sum_log_1, @sum_log_2, @sum_log_1 };
    params.zero = 0;
    params.constant = 1;
    params.reduce = true;
    params.reverse = false;
    params.name = 'sum_log';
    params.dimarg = 2;
end

try
    y = cvx_reduce_op( params, varargin{:} );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = sum_log_c(x)
y = sum(builtin('log',x),1)

function y = sum_log_1(x)
y = sum(log(x),1);

function y = sum_log_2( x )
y = size(x,1) * log(geo_mean(x,1));	

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
