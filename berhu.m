function y = berhu( x, M, t )

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
%       BERHU is convex and nonmonotonic in X and T; therefore, when used in 
%       CVX specifications, X and T must be affine. T must be real.

error( nargchk( 1, 3, nargin ) );
if nargin < 2,
    M = 1;
elseif ~isnumeric( M ),
    error( 'Second argument must be numeric.' );
elseif ~isreal( M ) | any( M( : ) <= 0 ),
    error( 'Second argument must be real and positive.' );
end

if nargin < 3,
    t = 1;
elseif ~isreal( t ),
    error( 'Third argument must be real.' );
elseif cvx_isconstant( t ) & nnz( cvx_constant( t ) <= 0 ),
    error( 'Third argument must be real and positive.' );
end

if ~cvx_isaffine( x ) | ~cvx_isaffine( t ),
    error( sprintf( 'Disciplined convex programming error:\n    HUBER is convex and nonmonotonic; its arguments must therefore be affine.' ) );
end

%
% Check sizes
%

sx = size( x ); xs = all( sx == 1 );
sM = size( M ); Ms = all( sM == 1 );
st = size( t ); ts = all( st == 1 );
if ~xs, sz = sx; elseif ~Ms, sz = sM; else sz = st; end
if ~( xs | isequal( sz, sx ) ) | ~( Ms | isequal( sz, sM ) ) | ~( ts | isequal( sz, st ) ),
   error( 'Sizes are incompatible.' );
end

%
% Compute result
%

y = abs( x ./ max(t,realmin) );
z = min( y, M );
y = t .* ( y + ( y - z ).^2 / (2*M) );
if nnz( t <= 0 ),
    if ts, 
        y = Inf * ones( size( y ) );
    else
        y( t <= 0 ) = Inf;
    end
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
