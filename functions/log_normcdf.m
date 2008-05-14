function y = log_normcdf( x )

%LOG_NORMCDF   Logarithm of the cumulative normal distribution.
%   Y = LOG_NORMCDF(X) is the logarithm of the CDF of the normal
%   distribution at the point X.
%
%                           1    / x
%       LOG_NORMCDF(X) = ------- |   exp(-t^2/2) dt
%                        sqrt(2) / -Inf
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

error(nargchk(1,1,nargin));
if ~isreal( x ),
    error( 'Argument must be real.' );
end
y = log(0.5*erfc(-x*sqrt(0.5)));

a =sqrt( [ 0.018102332171520
           0.011338501342044
           0.072727608432177
           0.184816581789135
           0.189354610912339
           0.023660365352785 ] );
b = [3 2.5 2 1 -1 -2]';
nb = length(b);


sx = size(x);
nx = prod(sx);
x  = reshape( x, 1, nx );
y  = b( :, ones(nx,1) ) - x( ones(nb,1), : );
y  = sparse(diag(a)) * y;
y  = - reshape( sum_square_pos( y ), sx );
