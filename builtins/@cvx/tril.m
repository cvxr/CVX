function x = tril( x, k )

%   Disciplined convex/geometric programming information for TRIL:
%       TRIL imposes no convexity restrictions on its arguments.

%
% Check inputs
%

if nargin < 2, k = 0; end
s = x.size_;
if length( s ) > 2,
    error( 'The first argument must be 2-D.' );
elseif ~isnumeric( k ) || length( k ) ~= 1,
    error( 'The second argument must be an integer scalar.' );
end

%
% Zero out the elements outside of the desired triangle
%

b = x.basis_;
b( :, ~tril(ones(s),k) ) = 0;
x = cvx( s, b );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
