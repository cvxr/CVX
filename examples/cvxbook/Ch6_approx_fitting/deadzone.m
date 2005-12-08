function y = deadzone( x, M )
error( nargchk( 1, 2, nargin ) );

% DEADZONE   Deadzone penalty function.
%     For a real or complex scalar X, DEADZONE(X) is the deadzone penalty
%     function applied to X: that is,
%
%         DEADZONE(X) = 0    if |X|<=1,
%                       |X|-1  if |X|>=1. 
%
%     DEADZONE(X,M) is the deadzone penalty function of halfwidth M applied to X;
%     that is, DEADZONE(X,M)=DEADZONE(X./M).
%
%     For matrices and N-D arrays, the penalty function is applied to each
%     element of X independently. M and X must be compatible in the same
%     sense as .*: one must be a scalar, or they must have identical size.
%
%     Disciplined convex programming information:
%         DEADZONE is convex and nonmonotonic; therefore, when used in CVX
%         specifications, its argument must be affine.

%
% Check types
%

if nargin < 2,
    M = 1;
elseif ~isnumeric( M ),
    error( 'Second argument must be numeric.' );
elseif ~isreal( M ) || any( M( : ) <= 0 ),
    error( 'Second argument must be real and positive.' );
end

%
% Check sizes
%

sx = size( x );
sM = size( M );
if all( sx == 1 ),
   sz = sM;
elseif all( sM == 1 ) || isequal( sx, sM ),
   sz = sx;
else
   error( 'Sizes are incompatible.' );
end

%
% Compute result
%

y = max(abs(x) - 1, 0);

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
