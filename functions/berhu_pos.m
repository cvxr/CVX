function y = berhu_pos( x, M, t )

%BERHU_POS   Monotonic reverse Huber penalty function.
%   BERHU_POS(X) computes the monotonic version of the reverse Huber penalty
%
%       BERHU(X) = MAX(X,0)     if X<=1,
%                  (X.^2+1)./2  if X>=1.
%
%   BERHU_POS(X,M) computes the reverse Huber penalty with halfwidth M,
%   M.*BERHU_POS(X./M). M must be real and positive.
%
%   BERHU_POS(X,M,T) computes the reverse Huber penalty with halfwidth M
%   and concomitant scale T: 
%
%       BERHU__)S(X,M,T) = (M.*T).*BERHU_POS(X./(M.*T)) if T>0
%                          +Inf                         if T<=0
%
%   This form supports the joint estimation of regression coefficients and
%   scaling; c.f. Art B. Owen, "A robust hybrid of lasso and ridge regression",
%   techincal report, Department of Statistics, Stanford University, 2006: 
%       http://www-stat.stanford.edu/~owen/reports/hhu.pdf
%
%   For matrices and N-D arrays, the penalty function is applied to each
%   element of X independently. M and T must be compatible with X in the same
%   sense as .*: one must be a scalar, or they must have identical size.
%
%   Disciplined convex programming information:
%       BERHU_POS is jointly convex in X and T, nondecreasing in X, and
%       nonincreasing in T. Therefore, when used in CVX specifications, X
%       must be convex and T must be concave. T must be real.

%
% Check arguments
%

error( nargchk( 1, 3, nargin ) ); %#ok
if nargin < 2,
    M = 1;
elseif ~isreal( M ) || any( M( : ) <= 0 ),
    error( 'Second argument must be real and positive.' );
end
if nargin < 3,
    t = 1;
elseif ~isreal( t ),
    error( 'Third argument must be real.' );
end
sz = cvx_size_check( x, M, t );
if isempty( sz ),
    error( 'Sizes are incompatible.' );
end

%
% Compute result
%

y = max( x, 0 ) ./ max(t,realmin);
z = min( y, M );
y = t .* ( y + ( y - z ).^2 ./ (2*M) );
q = t <= 0;
if nnz( q ),
    if length(t) == 1, 
        y = Inf * ones( sz );
    else
        y( q ) = Inf;
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
