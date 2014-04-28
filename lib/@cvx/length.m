function n = length( x )
s = x.size_;
if all( s ),
   n = max( s );
else
   n = 0;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
