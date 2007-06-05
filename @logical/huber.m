function y = huber( x, M )
error( nargchk( 1, 2, nargin ) );

% HUBER   Huber penalty function.
%     For a real or complex scalar X, HUBER(X) is the Huber penalty
%     function applied to X: that is,
%
%         HUBER(X) = |X|^2   if |X|<=1,
%                    2|X|-1  if |X|>=1.
%
%     HUBER(X,M) is the Huber penalty function of halfwidth M applied to X;
%     that is, HUBER(X,M)=M.^2.*HUBER(X./M).
%
%     For matrices and N-D arrays, the penalty function is applied to each
%     element of X independently. M and X must be compatible in the same
%     sense as .*: one must be a scalar, or they must have identical size.
%
%     Disciplined convex programming information:
%         HUBER is convex and nonmonotonic; therefore, when used in CVX
%         specifications, its argument must be affine.

%
% Check types
%

if nargin < 2,
    M = 1;
elseif ~isnumeric( M ),
    error( 'Second argument must be numeric.' );
elseif ~isreal( M ) | any( M( : ) <= 0 ),
    error( 'Second argument must be real and positive.' );
end

%
% Check sizes
%

sx = size( x );
sM = size( M );
if all( sx == 1 ),
   sz = sM;
elseif all( sM == 1 ) | isequal( sx, sM ),
   sz = sx;
else
   error( 'Sizes are incompatible.' );
end

%
% Compute result
%

y = abs( x );
z = min( y, M );
y = z .* ( 2 * y - z );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
