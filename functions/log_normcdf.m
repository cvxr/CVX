function y = log_normcdf( x, approx )

%LOG_NORMCDF   Logarithm of the cumulative normal distribution.
%   Y = LOG_NORMCDF(X) is the logarithm of the CDF of the normal
%   distribution at the point X.
%
%                                1    / x
%       LOG_NORMCDF(X) = LOG( ------- |   exp(-t^2/2) dt )
%                             sqrt(2) / -Inf
%
%   For numeric X, LOG_NORMCDF(X) is computed using the equivalent 
%   expression LOG(0.5*ERFC(-X*SQRT(0.5))). When X is a CVX variable, a 
%   a piecewise quadratic *approximation* is employed instead. This
%   approximation gives good results when -4 <= x <= 4, and will be
%   improved in future releases of CVX.
%
%   For array values of X, the LOG_NORMCDF returns an array of identical
%   size with the calculation applied independently to each element.
%
%   X must be real.
%
%   Disciplined convex programming information:
%       LOG_NORMCDF is concave and nondecreasing in X. Therefore, when used
%       in CVX specifications, X must be concave.

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'real' }, { 'concave' } );
    P.funcs = { @log_normcdf_cnst, @log_normcdf_cncv };
end
if nargin == 2 && approx,
    P.funcs{1} = @log_normcdf_cncv;
else
    P.funcs{1} = @log_normcdf_cnst;
end
y = cvx_unary_op( P, x );

function y = log_normcdf_cnst( x )
y = log(0.5*erfc(-x*sqrt(0.5)));

function y = log_normcdf_cncv( x )
persistent a b
if isempty( a ),
    a = sqrt( [ 0.018102332171520
                0.011338501342044
                0.072727608432177
                0.184816581789135
                0.189354610912339
                0.023660365352785 ]' );
    b = [3 2.5 2 1 -1 -2];
end
y = bsxfun( @times, a, bsxfun( @minus, b, x ) );
y = - sum_square( pos( y ), 2 );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
