function y = isreal( x, full )

%Disciplined convex/geometric programming information for ISREAL:
%   ISREAL(X) may be freely applied to any CVX expression. It will
%   return TRUE for all real affine, convex, concave, log-convex,
%   and log-concave expressions, and FALSE for complex affine 
%   expressions.

y = x.basis_;
if isreal( y ),
    if nargin > 1 && full,
        y = true( x.size_ );
    else
        y = true;
    end
elseif nargin > 1 && full,
    y = ~any( imag( y ), 1 );
else
    y = nnz( imag( y ) ) == 0;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
