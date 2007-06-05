function y = berhu( x, M, t )
error( nargchk( 1, 3, nargin ) );

% BERHU   Reverse Huber penalty function.
%     For a real or complex scalar X, BERHU(X) is the reverse Huber penalty
%     function applied to X: that is,
%
%         BERHU(X) = |X|          if |X|<=1,
%                    (|X|^2+1)/2  if |X|>=1.
%
%     BERHU(X,M) is the Huber penalty function of halfwidth M applied to X;
%     that is, BERHU(X,M)=M.*BERHU(X./M).
%
%     BERHU(X,M,T) computes T.*BERHU(X./T,M), the perspective transformation
%     of BERHU(X,M). This is useful for solving regression problems with
%     concomitant scale. T is constrained to be nonnegative.
%
%     For matrices and N-D arrays, the penalty function is applied to each
%     element of X independently. M and X must be compatible in the same
%     sense as .*: one must be a scalar, or they must have identical size.
%
%     Disciplined convex programming information:
%         BERHU is convex and nonmonotonic; therefore, when used in CVX
%         specifications, its argument must be affine.

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

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
