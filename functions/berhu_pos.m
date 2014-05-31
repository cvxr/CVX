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

persistent P
if isempty( P ),
    P.map = cvx_remap( ...
        { { 'any' }, { 'nonpositive' } }, ...
        { { 'real' } }, ...
        { { 'convex' }, { 'concave' } }, [0,1,2] );
    P.funcs = { @berhu_pos_c, @berhu_pos_nc };
    P.constant = 1;
    P.name = 'berhu_pos';
end
if nargin < 2,
    M = 1;
elseif ~( isnumeric(M) && numel(M)==1 && isreal(M) && M>0 ),
    cvx_throw( Second argument must be a positive scalar.' );
end
if nargin < 3,
    t = 1;
end
y = cvx_binary_op( P, x, t, M );

function y = berhu_pos_c( x, t, M )
y = max( x, 0 ) ./ max(t,realmin);
z = min( y, M );
y = t .* ( y + ( y - z ).^2 ./ (2*M) );

function cvx_optval = berhu_pos_nc( x, t, M ) %#ok
cvx_begin
    variable v( sz ) nonnegative
    variable w( sz ) nonnegative
    minimize( quad_over_lin( w, t, 0 ) ./ ( 2 * M ) + w + v )
    x <= w + v; %#ok
    v <= M .* t; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
