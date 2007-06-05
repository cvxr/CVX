function y = isreal( x, full )
y = x.basis_;
if nargin > 1 & full,
    y = any( imag( y ), 1 );
else
    y = isreal( x.basis_ );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
