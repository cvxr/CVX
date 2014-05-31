function x = transpose( x )

%   Disciplined convex/geometric programming information for TRANSPOSE:
%      The transpose operation may be applied to CVX variables without
%      restriction.

s = x.size_;
if length( s ) > 2,
    cvx_throw( 'Transpose of an ND array is not defined.' );
elseif any( s > 1 ),
    b = x.basis_;
    if all( s > 1 ),
        b = b(:,reshape( 1 : prod( s ), s )'); 
    end
    x = cvx( [ s(2), s(1) ], b );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
