function y = berhu( x, varargin )

%BERHU   Reverse Huber penalty function.
%   BERHU(X) computes the reverse Huber penalty function
%
%       BERHU(X) = ABS(X)            if ABS(X)<=1,
%                  (ABS(X).^2+1)./2  if ABS(X)>=1.
%
%   BERHU(X,M) computes the reverse Huber penalty function with halfwidth M,
%   M.*BERHU(X./M). M must be real and positive.
%
%   BERHU(X,M,T) computes the reverse Huber penalty function with halfwidth M
%   and concomitant scale T: 
%
%       BERHU(X,M,T) = (M.*T).*BERHU(X./(M.*T))     if T>0
%                      +Inf                         if T<=0
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
%       BERHU is jointly convex in X and T, nondecreasing for positive X,
%       and nonincreasing in T. Therefore, when used in CVX specifications, 
%       X must be affine or nonnegative convex, and T must be concave or
%       affine. T must be real.

if ~cvx_isaffnnc( x ),
    cvx_throw( 'Disciplined convex programming error:\n    X must be affine or nonnegative convex.' );
end
y = berhu_pos( abs( x ), varargin{:} );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
